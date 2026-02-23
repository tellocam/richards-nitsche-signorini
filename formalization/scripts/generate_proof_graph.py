#!/usr/bin/env python3
"""
generate_proof_graph.py

Generates the proof-chain dependency graph for RichardsNitscheSignorini/Main.lean
and renders it as a PDF using graphviz.

Usage:
    python generate_proof_graph.py [--output PATH] [--dot-only]

Output:  ../graphs/proof_chain.pdf   (default)
         ../graphs/proof_chain.dot   (always written)
Requires: graphviz installed (dot command in PATH) for PDF rendering

The graph is hardcoded from the known proof structure in Main.lean.
Node categories:
  - constitutive axioms  → red     (#ff6b6b)  / dashed border
  - established results  → blue    (#4ecdc4)  / solid
  - novel results        → green   (#45b7d1)  / solid
  - proved theorems      → gold    (#f9ca24)  / solid
  - main theorems        → dark green (#2ecc71) / bold
"""

import argparse
import os
import shutil
import subprocess
import sys
import textwrap


# ---------------------------------------------------------------------------
# Graph definition
# ---------------------------------------------------------------------------

def build_dot() -> str:
    """Return a DOT string encoding the Main.lean proof-chain graph.

    Layout: an invisible separator node + staggered invisible edges force
    the two disconnected components (error analysis, constitutive properties)
    to stack vertically instead of side-by-side.
    """

    dot = textwrap.dedent("""\
        digraph proof_chain {
            // ----------------------------------------------------------------
            // Graph — compact stacked layout via invisible separator
            // ----------------------------------------------------------------
            graph [
                label = "RichardsNitscheSignorini — Proof Chain\\nMain.lean → Theorem 4.10"
                labelloc = "t"
                fontsize = 14
                fontname = "Helvetica"
                rankdir = "TB"
                bgcolor = "white"
                splines = "true"
                nodesep = 0.25
                ranksep = 0.5
                pad = 0.3
                newrank = true
            ]

            node [fontname = "Helvetica" fontsize = 8 style = "filled"
                  shape = "box" margin = "0.06,0.04" height = 0.01 width = 0.01]
            edge [color = "#555555" arrowsize = 0.5]

            // ================================================================
            // LEGEND — single HTML-table node
            // ================================================================
            legend [shape = "plaintext" label = <
                <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="4" CELLPADDING="3">
                <TR>
                  <TD BGCOLOR="#ff6b6b" STYLE="dashed" COLOR="#cc0000"><FONT FACE="Helvetica" POINT-SIZE="8">Constitutive axiom</FONT></TD>
                  <TD BGCOLOR="#4ecdc4" COLOR="#2a9d8f"><FONT FACE="Helvetica" POINT-SIZE="8">Established result</FONT></TD>
                  <TD BGCOLOR="#45b7d1" COLOR="#1a7fa0"><FONT FACE="Helvetica" POINT-SIZE="8">Novel result</FONT></TD>
                  <TD BGCOLOR="#f9ca24" COLOR="#c09000"><FONT FACE="Helvetica" POINT-SIZE="8">Proved theorem</FONT></TD>
                  <TD BGCOLOR="#2ecc71" COLOR="#1a8a4a"><FONT FACE="Helvetica" POINT-SIZE="8" COLOR="white">Main theorem</FONT></TD>
                </TR>
                </TABLE>
            >]

            // ================================================================
            // ERROR ANALYSIS COMPONENT (top half)
            // ================================================================

            // -- Error-analysis established results (teal) --
            prop_4_7_splitting_error   [label = "Prop 4.7\\nsplitting O(√(Δt/Gr))"  fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_8_godunov_sub       [label = "Prop 4.8\\nGodunov sub O(h^½)"     fillcolor = "#4ecdc4" color = "#2a9d8f"]
            diffusion_error_Oh         [label = "Diffusion\\nerror O(h)"             fillcolor = "#4ecdc4" color = "#2a9d8f"]

            // -- Novel result (blue) --
            prop_4_9_godunov_sat [label = "Prop 4.9 [Novel]\\nGodunov sat O(h^(m/(m+1)))" fillcolor = "#45b7d1" color = "#1a7fa0"]

            // -- Proved intermediate theorems (gold) --
            splitting_error_vanishes   [label = "Split vanishes\\nerror → 0"      fillcolor = "#f9ca24" color = "#c09000" penwidth = 2]
            splitting_improves_with_Gr [label = "Split ↓ w/ Gr\\n↑Gr → ↓error"   fillcolor = "#f9ca24" color = "#c09000" penwidth = 2]
            sat_rate_worse_than_sub    [label = "Sat rate worse\\nm/(m+1) < 1/2"  fillcolor = "#f9ca24" color = "#c09000" penwidth = 2]

            // -- Main theorems (green, bold) --
            thm_4_10_combined_sub [
                label = "Thm 4.10(a)\\nsub-saturated\\nO(min(Δt,√(Δt/Gr)) + h^½ + h)"
                fillcolor = "#2ecc71" fontcolor = "white" color = "#1a8a4a"
                style = "filled,bold" penwidth = 3 fontsize = 9
            ]
            thm_4_10_combined_sat [
                label = "Thm 4.10(b)\\nsaturated\\nO(min(Δt,√(Δt/Gr)) + h^(m/(m+1)) + h)"
                fillcolor = "#2ecc71" fontcolor = "white" color = "#1a8a4a"
                style = "filled,bold" penwidth = 3 fontsize = 9
            ]
            spatial_rate_sublinear [
                label = "Sublinear rate\\n0 < m < 1"
                fillcolor = "#2ecc71" fontcolor = "white" color = "#1a8a4a"
                style = "filled,bold" penwidth = 3 fontsize = 9
            ]
            splitting_accuracy_improves_with_degeneracy [
                label = "Split accuracy\\n↑Gr → ↑accuracy"
                fillcolor = "#2ecc71" fontcolor = "white" color = "#1a8a4a"
                style = "filled,bold" penwidth = 3 fontsize = 9
            ]

            // ================================================================
            // Invisible separator between the two components
            // ================================================================
            separator [style = invis width = 0 height = 0 fixedsize = true label = ""]

            // ================================================================
            // CONSTITUTIVE PROPERTIES COMPONENT (bottom half)
            // ================================================================

            // -- Constitutive axioms (red, dashed) --
            K_Se_zero              [label = "K(0) = 0\\nMualem 1976"       fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_one               [label = "K(1) = Ks\\nMualem 1976"      fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_strictMono        [label = "K strict mono\\nvG 1980"      fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_continuousOn      [label = "K continuous\\nvG 1980"       fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_nonneg            [label = "K ≥ 0\\nMualem 1976"          fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_holder_at_one     [label = "Hölder at sat\\neq. (4.20)"   fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_asymp_dry         [label = "K asymp dry\\neq. (4.11)"     fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            C_Se_asymp_dry         [label = "C asymp dry\\neq. (4.13)"     fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_deriv_zero_at_dry [label = "K' → 0 dry\\neq. (4.12)"     fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]
            K_Se_deriv_infty_at_sat[label = "K' → ∞ sat\\neq. (4.21)"     fillcolor = "#ff6b6b" color = "#cc0000" style = "filled,dashed"]

            // -- Constitutive established results (teal) --
            prop_4_1_dry_end_asymp     [label = "Prop 4.1\\nD ~ Se^β dry"           fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_1_D_tends_zero      [label = "Cor 4.1\\nD → 0 dry"              fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_2_saturation_blowup [label = "Prop 4.2\\nD blow-up sat"          fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_2_D_tends_top       [label = "Cor 4.2\\nD → ∞ sat"              fillcolor = "#4ecdc4" color = "#2a9d8f"]
            cor_4_6a_Kinv_blowup       [label = "Cor 4.6a\\nK⁻¹ blow-up"           fillcolor = "#4ecdc4" color = "#2a9d8f"]
            cor_4_6b_condition_number   [label = "Cor 4.6b\\ncond. number"          fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_4a_existence        [label = "Prop 4.4a\\nexistence"              fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_4b_L1_contraction   [label = "Prop 4.4b\\nL1 contraction"        fillcolor = "#4ecdc4" color = "#2a9d8f"]
            prop_4_5_front_ratio       [label = "Prop 4.5\\nfront ratio"             fillcolor = "#4ecdc4" color = "#2a9d8f"]
            K_Se_lipschitz_on_compact  [label = "K Lipschitz\\ncompact [0,1)"        fillcolor = "#4ecdc4" color = "#2a9d8f"]

            // -- Constitutive proved theorems (gold) --
            prop_4_3_RH_speed          [label = "Prop 4.3\\nRH speed s = Ks/Δθ"  fillcolor = "#f9ca24" color = "#c09000" penwidth = 2]
            prop_4_3_lax_entropy       [label = "Prop 4.3\\nLax entropy"          fillcolor = "#f9ca24" color = "#c09000" penwidth = 2]
            prop_4_4c_uniqueness       [label = "Prop 4.4c\\nuniqueness"          fillcolor = "#4ecdc4" color = "#2a9d8f"]
            K_inv_at_one               [label = "K⁻¹(1) = 1/Ks\\nfrom K_Se_one"  fillcolor = "#f9ca24" color = "#c09000" penwidth = 2]

            // ================================================================
            // EDGES
            // ================================================================

            // Axioms → Established / Proved
            K_Se_asymp_dry     -> prop_4_1_dry_end_asymp
            C_Se_asymp_dry     -> prop_4_1_dry_end_asymp
            prop_4_1_dry_end_asymp -> prop_4_1_D_tends_zero
            K_Se_holder_at_one -> cor_4_6a_Kinv_blowup
            K_Se_holder_at_one -> cor_4_6b_condition_number
            prop_4_2_saturation_blowup -> prop_4_2_D_tends_top
            K_Se_zero -> prop_4_3_RH_speed
            K_Se_one  -> prop_4_3_RH_speed
            K_Se_one  -> K_inv_at_one

            // Proved chains
            prop_4_3_RH_speed      -> prop_4_3_lax_entropy
            prop_4_4b_L1_contraction -> prop_4_4c_uniqueness

            // Splitting chain
            prop_4_7_splitting_error -> splitting_error_vanishes
            prop_4_7_splitting_error -> splitting_improves_with_Gr

            // Godunov rate comparison
            prop_4_8_godunov_sub -> sat_rate_worse_than_sub
            prop_4_9_godunov_sat -> sat_rate_worse_than_sub

            // → Main theorem 4.10(a)
            prop_4_7_splitting_error -> thm_4_10_combined_sub
            prop_4_8_godunov_sub     -> thm_4_10_combined_sub
            diffusion_error_Oh       -> thm_4_10_combined_sub

            // → Main theorem 4.10(b)
            prop_4_7_splitting_error -> thm_4_10_combined_sat
            prop_4_9_godunov_sat     -> thm_4_10_combined_sat
            diffusion_error_Oh       -> thm_4_10_combined_sat

            // → Assembly corollaries
            splitting_improves_with_Gr -> splitting_accuracy_improves_with_degeneracy

            // ================================================================
            // LAYOUT: Vertical stacking via separator
            // ================================================================

            // Pin legend to top
            { rank = min; legend }

            // Pin error-analysis main theorems together
            { rank = same; thm_4_10_combined_sub; thm_4_10_combined_sat;
                           spatial_rate_sublinear; splitting_accuracy_improves_with_degeneracy }

            // Invisible chain: error bottom → separator → constitutive top
            thm_4_10_combined_sub -> separator [style = invis minlen = 1]
            thm_4_10_combined_sat -> separator [style = invis minlen = 1]
            spatial_rate_sublinear -> separator [style = invis minlen = 1]
            splitting_accuracy_improves_with_degeneracy -> separator [style = invis minlen = 1]

            // Separator fans out to constitutive root nodes — staggered across 3 rows
            // Row A (4 nodes): axioms with children + saturation blowup
            separator -> K_Se_zero              [style = invis minlen = 1]
            separator -> K_Se_one               [style = invis minlen = 1]
            separator -> K_Se_asymp_dry         [style = invis minlen = 1]
            separator -> C_Se_asymp_dry         [style = invis minlen = 1]

            // Row B (5 nodes): pushed 1 rank below row A
            K_Se_zero  -> K_Se_holder_at_one          [style = invis]
            K_Se_one   -> prop_4_2_saturation_blowup  [style = invis]
            K_Se_asymp_dry  -> prop_4_4b_L1_contraction [style = invis]
            C_Se_asymp_dry  -> K_Se_strictMono        [style = invis]

            // Row C (6 nodes): pushed 1 rank below row B
            K_Se_holder_at_one -> K_Se_continuousOn   [style = invis]
            prop_4_2_saturation_blowup -> K_Se_nonneg [style = invis]
            prop_4_4b_L1_contraction -> K_Se_deriv_zero_at_dry [style = invis]
            K_Se_strictMono -> K_Se_deriv_infty_at_sat [style = invis]

            // Row D (remaining): pushed 1 rank below row C
            K_Se_continuousOn -> prop_4_4a_existence  [style = invis]
            K_Se_nonneg -> prop_4_5_front_ratio       [style = invis]
            K_Se_deriv_zero_at_dry -> K_Se_lipschitz_on_compact [style = invis]
        }
    """)
    return dot


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate proof-chain dependency graph for RichardsNitscheSignorini/Main.lean"
    )
    parser.add_argument(
        "--output",
        default=None,
        help="Output PDF path (default: ../graphs/proof_chain.pdf relative to this script)"
    )
    parser.add_argument(
        "--dot-only",
        action="store_true",
        help="Write the DOT file and exit without attempting PDF rendering"
    )
    args = parser.parse_args()

    # Determine paths
    script_dir = os.path.dirname(os.path.abspath(__file__))
    graphs_dir = os.path.join(os.path.dirname(script_dir), "graphs")
    os.makedirs(graphs_dir, exist_ok=True)

    dot_path = os.path.join(graphs_dir, "proof_chain.dot")
    pdf_path = args.output if args.output else os.path.join(graphs_dir, "proof_chain.pdf")

    # Write DOT file (always)
    dot_content = build_dot()
    with open(dot_path, "w", encoding="utf-8") as f:
        f.write(dot_content)
    print(f"DOT file written to: {dot_path}")

    if args.dot_only:
        print("--dot-only specified; skipping PDF rendering.")
        return

    # Try to render PDF with graphviz (optional — prints install hint if absent)
    if shutil.which("dot") is None:
        print("NOTE: graphviz 'dot' not found in PATH; DOT file written but PDF not rendered.",
              file=sys.stderr)
        print("Install graphviz to render PDF:", file=sys.stderr)
        print("  Fedora/RHEL:   sudo dnf install graphviz", file=sys.stderr)
        print("  Ubuntu/Debian: sudo apt install graphviz", file=sys.stderr)
        print("  macOS:         brew install graphviz", file=sys.stderr)
        sys.exit(0)

    # Render directly with dot — the DOT file already contains invisible
    # separator edges and stagger chains that produce a compact layout.
    # (unflatten is not used: it widens the stacked layout.)
    print(f"Rendering {dot_path} -> {pdf_path} ...")
    result = subprocess.run(
        ["dot", "-Tpdf", "-o", pdf_path, dot_path],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print("ERROR: graphviz dot failed:", file=sys.stderr)
        print(result.stderr, file=sys.stderr)
        sys.exit(result.returncode)

    print(f"\nDone. Proof chain graph written to: {pdf_path}")


if __name__ == "__main__":
    main()
