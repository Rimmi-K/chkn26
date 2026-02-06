import re
import numpy as np
import pandas as pd
from cobra import Model
from typing import List, Dict, Tuple, Union, Optional


def build_gpr_rule_map(model: Model) -> Dict[str, str]:
    """Map reaction ID to GPR string"""
    return {rxn.id: (rxn.gene_reaction_rule or "") for rxn in model.reactions}


def build_reaction_pathways_map(model: Model) -> Dict[str, str]:
    """Map reaction ID to subsystem/pathway"""
    return {rxn.id: (rxn.subsystem or "Unknown") for rxn in model.reactions}


Token = str
AST = Union[str, Tuple[str, "AST", "AST"]]

_TOKEN_RE = re.compile(r'\s*(\(|\)|and|or|[^()\s]+)\s*', flags=re.IGNORECASE)


def _tokenize(rule: str) -> List[Token]:
    if not rule:
        return []
    toks: List[Token] = []
    for m in _TOKEN_RE.finditer(rule):
        t = m.group(1)
        tl = t.lower()
        if tl in ("and", "or"):
            toks.append(tl)
        elif t in ("(", ")"):
            toks.append(t)
        else:
            toks.append(t)
    return toks


def _parse_gpr_to_ast(tokens: List[Token]) -> AST:
    """
    Recursive descent parser for GPR grammar:
      expr := term ( 'or' term )*
      term := factor ( 'and' factor )*
      factor := GENE | '(' expr ')'
    """
    idx = 0

    def peek():
        nonlocal idx
        return tokens[idx] if idx < len(tokens) else None

    def consume(expected=None):
        nonlocal idx
        tok = tokens[idx] if idx < len(tokens) else None
        if expected is not None and tok != expected:
            raise ValueError(f"Expected {expected!r}, got {tok!r}")
        idx += 1
        return tok

    def parse_factor() -> AST:
        tok = peek()
        if tok == '(':
            consume('(')
            node = parse_expr()
            consume(')')
            return node
        if tok is None:
            raise ValueError("Unexpected end of input in GPR")
        consume()
        return tok

    def parse_term() -> AST:
        node = parse_factor()
        while True:
            tok = peek()
            if tok == 'and':
                consume('and')
                right = parse_factor()
                node = ('AND', node, right)
            else:
                break
        return node

    def parse_expr() -> AST:
        node = parse_term()
        while True:
            tok = peek()
            if tok == 'or':
                consume('or')
                right = parse_term()
                node = ('OR', node, right)
            else:
                break
        return node

    ast = parse_expr()
    if idx != len(tokens):
        raise ValueError("Trailing tokens in GPR rule")
    return ast


def _ast_to_dnf(node: AST) -> List[List[str]]:
    """Convert AST to Disjunctive Normal Form (list of AND-clauses)"""
    if isinstance(node, str):
        return [[node]]
    op, left, right = node
    if op == 'OR':
        return _ast_to_dnf(left) + _ast_to_dnf(right)
    if op == 'AND':
        L = _ast_to_dnf(left)
        R = _ast_to_dnf(right)
        out: List[List[str]] = []
        for a in L:
            for b in R:
                out.append(a + b)
        return out
    raise ValueError(f"Unknown op {op}")


def parse_gpr_to_dnf(rule: str) -> List[List[str]]:
    """
    Convert GPR string to DNF clauses.
    Example: '(G1 and (G2 or G3)) or G4' -> [['G1','G2'], ['G1','G3'], ['G4']]
    """
    rule = (rule or "").strip()
    if not rule:
        return []
    tokens = _tokenize(rule)
    ast = _parse_gpr_to_ast(tokens)
    dnf = _ast_to_dnf(ast)

    norm: List[List[str]] = []
    for clause in dnf:
        seen = set()
        cleaned: List[str] = []
        for g in clause:
            if g not in seen:
                seen.add(g)
                cleaned.append(g)
        if cleaned:
            norm.append(cleaned)
    return norm


def build_gpr_clauses(model: Model) -> Dict[str, List[List[str]]]:
    """Map reaction ID to DNF clauses (list of AND-clauses)"""
    return {rxn.id: parse_gpr_to_dnf(rxn.gene_reaction_rule or "") for rxn in model.reactions}


def reaction_log2fc_via_gpr(
    clauses: Dict[str, List[List[str]]],
    deg_df: pd.DataFrame,
    missing_strategy: str = "impute",
    impute_value: float = 0.0,
) -> Dict[str, float]:
    """
    Calculate reaction log2FC from GPR using FASTCORMICS logic:
      - Within AND-clause: minimum across genes (weakest link)
      - Between OR-clauses: maximum across isoforms (strongest isoform)

    Parameters
    ----------
    clauses : dict
        Reaction ID to DNF clauses mapping
    deg_df : pd.DataFrame
        DEG table with columns 'gene' and 'log2fc'
    missing_strategy : str
        "skip" (drop clause if gene missing) or "impute" (use impute_value)
    impute_value : float
        Value to use for missing genes when strategy is "impute"

    Returns
    -------
    dict
        Reaction ID to log2FC mapping
    """
    if "gene" not in deg_df.columns or "log2fc" not in deg_df.columns:
        raise ValueError("deg_df must contain columns: 'gene', 'log2fc'")

    lfc_map = dict(zip(deg_df["gene"], deg_df["log2fc"]))
    rxn_lfc: Dict[str, float] = {}

    for rid, dnf in clauses.items():
        clause_aggs: List[float] = []

        for clause in dnf:
            vals: List[float] = []
            missing_any = False

            for g in clause:
                val = lfc_map.get(g, None)
                if val is None or pd.isna(val):
                    if missing_strategy == "impute":
                        vals.append(float(impute_value))
                    elif missing_strategy == "skip":
                        missing_any = True
                        break
                    else:
                        raise ValueError(f"unknown missing_strategy: {missing_strategy}")
                else:
                    vals.append(float(val))

            if missing_any:
                continue
            if vals:
                clause_aggs.append(min(vals))

        if clause_aggs:
            rxn_lfc[rid] = float(max(clause_aggs))

    return rxn_lfc



def reaction_pvalues_via_gpr(
    clauses: Dict[str, List[List[str]]],
    deg_df: pd.DataFrame,
    p_col: str = "padj",
    require_all_genes: bool = True
) -> Dict[str, float]:

    if not {"gene", p_col}.issubset(deg_df.columns):
        raise ValueError(f"deg_df must contain columns: 'gene', '{p_col}'")

    p_map = dict(zip(deg_df["gene"], deg_df[p_col]))
    rxn_p: Dict[str, float] = {}

    for rid, dnf in clauses.items():
        if not dnf:
            continue

        clause_pvals: List[float] = []

        for clause in dnf:
            genes = [g for g in clause if g in p_map and pd.notna(p_map[g])]
            if require_all_genes and len(genes) != len(clause):
                continue
            if not genes:
                continue

            pvals = [float(p_map[g]) for g in genes]
            pc = float(np.max(pvals))

            if pd.notna(pc):
                clause_pvals.append(pc)

        if not clause_pvals:
            continue

        prxn = _harmonic_mean_p(clause_pvals)
        rxn_p[rid] = prxn

    return rxn_p


def make_gpr_rule_with_values(rule: str,
                              rxn_genes: List[str],
                              lfc_map: Dict[str, float],
                              na_token: str = "NA",
                              fmt: str = "{:.2f}") -> str:
    """
    Replace gene IDs in GPR string with their numeric values.
    Replacement is done safely using word boundaries.
    """
    s = rule if rule else ""
    for gid in sorted(rxn_genes, key=len, reverse=True):
        val = lfc_map.get(gid, np.nan)
        rep = na_token if (val is None or pd.isna(val)) else fmt.format(float(val))
        s = re.sub(rf'(?<![\w]){re.escape(gid)}(?![\w])', rep, s)
    return s
