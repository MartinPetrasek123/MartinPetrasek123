#!/usr/bin/env python3
"""Build the standalone RFG-R completion paper as a verified PDF."""

from __future__ import annotations

from pathlib import Path

from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm
from reportlab.platypus import (
    Image,
    KeepTogether,
    PageBreak,
    Paragraph,
    Preformatted,
    SimpleDocTemplate,
    Spacer,
    Table,
    TableStyle,
)


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "paper" / "R_Universe_RFG_R_Completion.pdf"
FIGURES = ROOT / "generated" / "figures"


def page_number(canvas, doc):
    canvas.saveState()
    canvas.setStrokeColor(colors.HexColor("#bcc6cc"))
    canvas.line(doc.leftMargin, 1.45 * cm, A4[0] - doc.rightMargin, 1.45 * cm)
    canvas.setFillColor(colors.HexColor("#46535a"))
    canvas.setFont("Helvetica", 8)
    canvas.drawString(doc.leftMargin, 0.92 * cm, "R-Universe RFG-R Completion")
    canvas.drawRightString(A4[0] - doc.rightMargin, 0.92 * cm, f"{doc.page}")
    canvas.restoreState()


def styles():
    base = getSampleStyleSheet()
    return {
        "title": ParagraphStyle(
            "TitleCustom",
            parent=base["Title"],
            fontName="Helvetica-Bold",
            fontSize=20,
            leading=25,
            alignment=TA_CENTER,
            spaceAfter=8,
            textColor=colors.HexColor("#142b3a"),
        ),
        "author": ParagraphStyle(
            "Author",
            parent=base["Normal"],
            fontName="Helvetica",
            fontSize=10.5,
            leading=14,
            alignment=TA_CENTER,
            spaceAfter=20,
            textColor=colors.HexColor("#46535a"),
        ),
        "h1": ParagraphStyle(
            "H1",
            parent=base["Heading1"],
            fontName="Helvetica-Bold",
            fontSize=14,
            leading=18,
            spaceBefore=13,
            spaceAfter=7,
            textColor=colors.HexColor("#0d4057"),
        ),
        "h2": ParagraphStyle(
            "H2",
            parent=base["Heading2"],
            fontName="Helvetica-Bold",
            fontSize=11,
            leading=14,
            spaceBefore=9,
            spaceAfter=5,
            textColor=colors.HexColor("#0d4057"),
        ),
        "body": ParagraphStyle(
            "Body",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=9.5,
            leading=13.6,
            alignment=TA_JUSTIFY,
            spaceAfter=7,
        ),
        "abstract": ParagraphStyle(
            "Abstract",
            parent=base["BodyText"],
            fontName="Helvetica-Oblique",
            fontSize=9.5,
            leading=13.6,
            alignment=TA_JUSTIFY,
            leftIndent=0.55 * cm,
            rightIndent=0.55 * cm,
            spaceAfter=10,
        ),
        "eq": ParagraphStyle(
            "Equation",
            parent=base["Code"],
            fontName="Courier",
            fontSize=8.8,
            leading=12.1,
            leftIndent=0.7 * cm,
            rightIndent=0.4 * cm,
            backColor=colors.HexColor("#f4f7f8"),
            borderColor=colors.HexColor("#d8e1e5"),
            borderWidth=0.4,
            borderPadding=5,
            spaceBefore=3,
            spaceAfter=8,
        ),
        "caption": ParagraphStyle(
            "Caption",
            parent=base["BodyText"],
            fontName="Helvetica-Oblique",
            fontSize=8.3,
            leading=11,
            alignment=TA_CENTER,
            spaceAfter=10,
            textColor=colors.HexColor("#46535a"),
        ),
        "ref": ParagraphStyle(
            "Reference",
            parent=base["BodyText"],
            fontName="Helvetica",
            fontSize=8.5,
            leading=11.3,
            leftIndent=0.45 * cm,
            firstLineIndent=-0.45 * cm,
            spaceAfter=4,
        ),
    }


def P(text: str, style: ParagraphStyle):
    return Paragraph(text, style)


def E(text: str, style: ParagraphStyle):
    return Preformatted(text, style)


def image(path: Path, caption: str, st: dict[str, ParagraphStyle]):
    item = Image(str(path), width=16.1 * cm, height=10.3 * cm)
    item.hAlign = "CENTER"
    return KeepTogether([item, Spacer(1, 2), P(caption, st["caption"])])


def main() -> None:
    st = styles()
    doc = SimpleDocTemplate(
        str(OUT),
        pagesize=A4,
        leftMargin=2.05 * cm,
        rightMargin=2.05 * cm,
        topMargin=1.8 * cm,
        bottomMargin=2.0 * cm,
        title="R-Universe RFG-R Completion",
        author="Martin Petrasek",
    )
    story = []
    story += [
        P("R-Universe RFG-R", st["title"]),
        P("A regular multiscale completion with a matter, CMB, and PPN likelihood definition", st["title"]),
        P("Martin Petrasek", st["author"]),
        P("<b>Abstract.</b> RFG-R is a precisely specified low-energy completion of the R-Universe preferred-foliation branch. It regularizes the X=0 weak-field singularity while retaining the reconstructed cosmological background and luminal tensor propagation to O((epsilon/X)^p). The completion has two explicit domains: a cosmological relational action and a locally GR-matched weak-field EFT. The potential is reconstructed by quadrature from the prescribed branch rather than postulated. The local domain has gamma=beta=1 and alpha1=alpha2=0. The cosmological domain supplies the exact action derivatives needed by a full 3+1 Einstein-Boltzmann solver, with stability gating before every CMB likelihood evaluation. This paper defines the model and the falsification pipeline; it does not claim a data preference before the stated likelihood is executed.", st["abstract"]),
        P("<b>Keywords:</b> modified gravity; preferred foliation; CMB; post-Newtonian tests; effective field theory; standard sirens", st["body"]),
        P("1. The completion problem", st["h1"]),
        P("The R-Universe homogeneous branch is predictive for H(a), the effective relational equation of state, and the tensor distance ratio. Its original response Q(X)=1-A X^(-theta), however, is singular in the X=0 limit used by a standard static weak-field expansion. A genuine empirical theory therefore needs a definition of its local limit as well as a source for the linear matter and radiation equations. RFG-R supplies both definitions. It is deliberately an effective low-energy model: cosmological geometry activates the relational response, while resolved local gravity is matched to General Relativity.", st["body"]),
        P("2. Regular cosmological response", st["h1"]),
        P("Let p be even and larger than theta, with reference values theta=1.6, p=4, and epsilon=10^-8. Define nu=1+theta/p and A=Omega_R0/(1+theta).", st["body"]),
        E("R_epsilon(X) = Omega_R0 X^(p+2)/(X^p+epsilon^p)^nu\nQ_epsilon(X) = 1 - A X^p/(X^p+epsilon^p)^nu", st["eq"]),
        P("The covariant preferred-foliation action is the original RFG action after replacing Q and V by Q_epsilon and V_epsilon. Matter remains minimally and universally coupled to the Jordan metric.", st["body"]),
        E("S_cos = (M_Pl^2/2) int d^4x sqrt(-g)\n        { Q_epsilon(X)[R3 + K_mn K^mn - K^2] + 2 H0^2 V_epsilon(X) }\n        + S_m[g,Psi]", st["eq"]),
        P("The FLRW lapse equation is fixed to be", st["body"]),
        E("E^2(a) = Omega_m0 a^-3 + Omega_r0 a^-4 + R_epsilon(E)", st["eq"]),
        P("and determines the potential without a new function:", st["body"]),
        E("F_epsilon = X^2 - R_epsilon - X^2 Q_epsilon - X^3 dQ_epsilon/dX\nV_epsilon(X) = -3 X int_0^X ds F_epsilon(s)/s^2", st["eq"]),
        P("The integral is finite because F_epsilon/s^2 is O(X^p) at the origin. A term proportional to X is an ADM boundary representative and has no field-equation content.", st["body"]),
        image(FIGURES / "regularization_recovery.png", "Figure 1. The regularized response becomes the original response immediately above the epsilon transition. On the cosmological branch the analytic difference is O((epsilon/X)^4).", st),
        P("3. Exact limits", st["h1"]),
        P("At X much larger than epsilon the regularized branch recovers the original one. With the reference epsilon, the analytic relative correction is below 2.5e-32 for X>=0.8; double precision evaluates the equality at about 1e-15. At X=0 the action is regular.", st["body"]),
        E("Q_epsilon = 1 - A epsilon^(-p-theta) X^p + O(X^(2p))\nV_epsilon = -[3 Omega_R0(p-theta)/((1+theta)(p+1)epsilon^(p+theta))] X^(p+2)\n              + O(X^(2p+2))", st["eq"]),
        P("For p=4 there is no physical correction to the quadratic Einstein-Hilbert ADM operator. This removes the X=0 singularity without altering the background or tensor branch that motivated RFG.", st["body"]),
        image(FIGURES / "fractional_densities.png", "Figure 2. Fractional densities of the reference RFG-R branch. The regularization is absent at visible plotting precision on this cosmological domain.", st),
        P("4. Local GR matching and PPN", st["h1"]),
        P("A single cosmological EFT is not automatically a local Solar-System EFT. RFG-R makes the crossover a model axiom. Define W=sqrt(abs(C_abcd C^abcd))/(H0/c)^2 and use a C-infinity switching function s(W) that equals one for W<=10^8 and zero for W>=10^9. The complete low-energy action is S_eff=S_GR+int s(W)[L_cos-L_GR]. FLRW has W=0, while a Schwarzschild Sun gives W=5.756e22 at one AU.", st["body"]),
        P("Inside the local matched domain s=0, so the action and all its variations are exactly those of GR. The PPN prediction is therefore gamma=1, beta=1, alpha1=0, and alpha2=0. The Cassini factor used by the likelihood is -2 ln L=[(gamma-1-2.1e-5)/(2.3e-5)]^2, yielding 0.833648 for the GR prediction.", st["body"]),
        P("5. Matter and CMB likelihood", st["h1"]),
        Spacer(1, 7),
        P("No compressed CMB distance prior is used. The CMB interface begins with the action-level ADM derivatives. For L=Q_epsilon(X)[R3+K_ij K^ij-K^2]+2H0^2 V_epsilon(X), the principal inputs include L_R3=Q_epsilon and L_Kij_R3=Q_epsilon,X gamma^ij/(3H0). The full kinetic Hessian follows by differentiating L_Kij once more. The repository generates all background functions and derivatives in eft_coefficients.csv, avoiding any differentiation of a singular power law.", st["body"]),
        P("For each sampled parameter point, the full 3+1 Einstein-Boltzmann backend must solve the scalar constraints together with photon, baryon, CDM, neutrino, recombination, and lensing evolution. It rejects a point before likelihood evaluation if the scalar constraint matrix becomes singular or if a ghost, gradient instability, or Q_T<=0 occurs. Accepted points export TT, TE, EE, lensing, P(k,z), f sigma8(z), and the standard-siren distance.", st["body"]),
        PageBreak(),
        P("6. Data protocol", st["h1"]),
        Spacer(1, 7),
        P("The baseline log likelihood is ln L=ln L_Planck18+ln L_lensing+ln L_BAO+ln L_SN+ln L_RSD+ln L_PPN+ln L_GW. Planck is evaluated with the official low-ell temperature and polarization, high-ell TT/TE/EE, and lensing packages; the solver is a full modified-gravity Einstein-Boltzmann solver rather than a quasi-static growth approximation. The public H-EFTCAMB architecture is suitable because it supports full linear modified-gravity evolution, 3+1 EFT functions, and stability tests. Cobaya can sample the joint likelihood after the action derivative module is connected.", st["body"]),
        P("The prediction is now sharply falsifiable. The mathematics of the completion and the allowed output set are closed. Whether the model is preferred, allowed, or ruled out is an empirical result of this explicitly defined likelihood, not a conclusion from background algebra alone.", st["body"]),
        P("7. Numerical verification", st["h1"]),
    ]
    table_data = [
        ["check", "reference result"],
        ["high-X response recovery", "5.329e-15 numerical maximum"],
        ["high-X Q recovery", "2.600e-15 numerical maximum"],
        ["potential reconstruction residual", "8.145e-06 finite-difference maximum"],
        ["background density closure", "1.554e-13 maximum error"],
        ["minimum Q_T over a in [1e-8,1e2]", "0.6153847956"],
        ["Solar W at 1 AU", "5.756e22; local GR domain"],
    ]
    table = Table(table_data, colWidths=[7.4 * cm, 8.1 * cm], repeatRows=1)
    table.setStyle(TableStyle([
        ("BACKGROUND", (0, 0), (-1, 0), colors.HexColor("#0d4057")),
        ("TEXTCOLOR", (0, 0), (-1, 0), colors.white),
        ("FONTNAME", (0, 0), (-1, 0), "Helvetica-Bold"),
        ("FONTNAME", (0, 1), (-1, -1), "Helvetica"),
        ("FONTSIZE", (0, 0), (-1, -1), 8.5),
        ("LEADING", (0, 0), (-1, -1), 11),
        ("GRID", (0, 0), (-1, -1), 0.35, colors.HexColor("#ced9de")),
        ("VALIGN", (0, 0), (-1, -1), "MIDDLE"),
        ("ROWBACKGROUNDS", (0, 1), (-1, -1), [colors.white, colors.HexColor("#f4f7f8")]),
        ("LEFTPADDING", (0, 0), (-1, -1), 6),
        ("RIGHTPADDING", (0, 0), (-1, -1), 6),
        ("TOPPADDING", (0, 0), (-1, -1), 5),
        ("BOTTOMPADDING", (0, 0), (-1, -1), 5),
    ]))
    story += [table, Spacer(1, 10)]
    story += [
        P("The checks are generated by scripts/run_all.sh. They test the action-defined completion rather than a fitted surrogate.", st["body"]),
        P("References", st["h1"]),
        P("[1] M. Petrasek, Relational Capacity Dynamics and its Preferred-Foliation Completion (2026).", st["ref"]),
        P("[2] B. Hu, M. Raveri, N. Frusciante, and A. Silvestri, Effective Field Theory of Cosmic Acceleration: an implementation in CAMB, Phys. Rev. D 89, 103530 (2014), arXiv:1312.5742.", st["ref"]),
        P("[3] G. Ye et al., H-EFTCAMB: A Cobaya-Integrated, Python-Wrapped Extension of EFTCAMB for Covariant Horndeski Gravity, arXiv:2603.01662 (2026).", st["ref"]),
        P("[4] Planck Collaboration, Planck 2018 results. V. CMB power spectra and likelihoods, Astron. Astrophys. 641, A5 (2020), arXiv:1907.12875.", st["ref"]),
        P("[5] B. Bertotti, L. Iess, and P. Tortora, A test of general relativity using radio links with the Cassini spacecraft, Nature 425, 374-376 (2003).", st["ref"]),
    ]
    doc.build(story, onFirstPage=page_number, onLaterPages=page_number)
    print(OUT)


if __name__ == "__main__":
    main()
