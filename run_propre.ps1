# =====================================================================
#  run_propre.ps1 -- le balayage d'equivariance recalcule A NEUF, dans
#  un seul etat du code.
#
#  CE QUE CE RUN-CI PORTE EN PLUS
#  ------------------------------
#  Les verdicts des 5.36 a 5.38. Trois d'entre eux etaient deja dans le
#  code du balayage -- la levee des gardes `len(c) == 1` -- et ne
#  demandaient qu'un balayage ; le quatrieme est neuf : le lieu de base
#  des 5.37 / 5.38 est desormais interroge la ou le certificat de
#  surjectivite echoue, et une ligne dont le lieu de base est DEMONTRE
#  sur Y sort ECARTEE (`fibre = false`) au lieu d'indeterminee.
#
#  POURQUOI UN NOUVEAU DOSSIER, ET PAS scan_wilson5
#  ------------------------------------------------
#  Modifier `equivariance_f.py` change l'empreinte du code. Reprendre
#  dans scan_wilson5 y ferait cohabiter deux versions -- exactement ce
#  que le 5.35 a coute (4 049 candidats ecartes a tort), et la recette
#  `retirer_lots.py --verifier` refuserait le fichier, a juste titre.
#  scan_wilson5 n'est pas suivi par git : on n'y touche pas, on ne le
#  --reset pas, et une garde en tete de ce script le refuse.
#
#  RIEN N'EST DETRUIT. On ecrit dans scan_wilson6 ; scan_wilson5 reste
#  en place comme base de comparaison. Pas de --reset, donc pas d'etape
#  destructive du tout.
#
#  Duree : ~29 h 30 a 7 coeurs pour le balayage de reference, plus le
#  cout du lieu de base sur les strates (1,3) et (2,3) -- mesure par
#  `compter_strates.py` AVANT de lancer. Fractionnable a volonte : le
#  checkpoint est ecrit apres chaque lot.
#
#      .\run_propre.ps1                  # premiere fois (47 tests + ancres)
#      .\run_propre.ps1 -SansTests       # sessions suivantes
#      .\run_propre.ps1 -SansTests -ControleFinal   # tous les lots 'T' faits
#
#  -Dossier / -Source pour viser autre chose ; -ControleLieu N regle
#  l'echantillon du controle negatif du lieu de base (0 = aucun, et le
#  bilan le dira).
# =====================================================================
param(
    [switch]$SansTests,
    [switch]$ControleFinal,
    [int]$Coeurs = 7,
    [string]$Dossier = 'scan_wilson6',
    [string]$Source  = 'scan_wilson5',
    [int]$ControleLieu = 20
)

$ErrorActionPreference = 'Stop'
$dossier = $Dossier
$source  = $Source
$t0 = Get-Date

# ---------------------------------------------------------------------
#  GARDE : ne jamais ecrire dans le dossier SOURCE.
#
#  `scan_wilson5` n'est pas suivi par git (le .gitignore ecarte scan_*),
#  il a coute une trentaine d'heures, et le 5.35 a etabli qu'un balayage
#  refait n'est pas garanti identique. L'ecraser est irreversible en
#  pratique. Cette garde ne remplace pas la prudence ; elle rattrape la
#  faute de frappe.
# ---------------------------------------------------------------------
if ($dossier -eq $source) {
    Write-Host "  ARRET : le dossier de sortie est le dossier source ($source)." -ForegroundColor Red
    Write-Host "  Le balayage porte les verdicts des 5.36 a 5.38 : il ECRIT." -ForegroundColor Red
    exit 1
}
if ($dossier -eq 'scan_wilson5') {
    Write-Host "  ARRET : scan_wilson5 est le fichier qui fait foi et git ne le protege pas." -ForegroundColor Red
    Write-Host "  Ecrire ailleurs (-Dossier scan_wilson6)." -ForegroundColor Red
    exit 1
}

function Etape($n, $titre) {
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "  ETAPE $n -- $titre" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
}

# ---------------------------------------------------------------------
# 0. Le dossier et son entree. `equivariance_f.py` lit
#    <dossier>\results_equivariant.jsonl et rien d'autre.
# ---------------------------------------------------------------------
if (-not (Test-Path $dossier)) {
    New-Item -ItemType Directory -Path $dossier | Out-Null
    Write-Host "  $dossier cree" -ForegroundColor DarkGray
}
$entree = Join-Path $dossier 'results_equivariant.jsonl'
if (-not (Test-Path $entree)) {
    Copy-Item (Join-Path $source 'results_equivariant.jsonl') $entree
    Write-Host "  entree copiee depuis $source" -ForegroundColor DarkGray
}

# ---------------------------------------------------------------------
# 1. Le socle. 45 tests.
# ---------------------------------------------------------------------
if (-not $SansTests) {
    Etape 1 "Tests de non-regression (47 tests, ~5 min)"
    python tests_regression.py
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  Tests en echec -- NE RIEN LANCER." -ForegroundColor Red
        exit 1
    }
    # Les ancres : le branchement reproduit-il les verdicts deja etablis ?
    # Sur les REFERENCES seules ici -- la comparaison avec le balayage vient
    # a l'etape 3, une fois qu'il y a quelque chose a comparer.
    python -u ancres_port.py
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  Les ancres des 5.36 a 5.38 ne tombent pas juste -- NE RIEN LANCER." -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "  (tests sautes : -SansTests)" -ForegroundColor DarkGray
}

# ---------------------------------------------------------------------
# 1 bis. Le controle d'orbite ne vaut que si le membre de controle et son
#        representant ont ete evalues contre le MEME fichier. Fractionner
#        le run casse cela en silence. A passer une fois tous les lots
#        'T' termines : les controles se rejouent contre un JSONL complet.
# ---------------------------------------------------------------------
if ($ControleFinal) {
    Etape "1 bis" "Rejeu des controles d'orbite"
    python -u retirer_lots.py $dossier --refaire-controles --appliquer
    if ($LASTEXITCODE -ne 0) { exit 1 }
}

# ---------------------------------------------------------------------
# 2. LA ligne, invariante d'une session a l'autre. Ne pas la modifier :
#    --cicy, --replier-orbites, --controle-orbites et --taille-lot entrent
#    dans l'empreinte du checkpoint, et une empreinte qui ne correspond
#    plus fait EFFACER le JSONL.
# ---------------------------------------------------------------------
Etape 2 "Equivariance de f, balayage complet (Ctrl-C quand vous voulez)"
try {
    python -u equivariance_f.py cicyquotients.m cicylist.txt $dossier `
        --replier-orbites -j $Coeurs --controle-lieu-de-base $ControleLieu |
        Tee-Object -Append -FilePath "${dossier}_equiv_f.log"
} finally {
    # -----------------------------------------------------------------
    # 3. Recette, executee meme apres un Ctrl-C. Sur un fichier propre,
    #    les trois doivent tomber :
    #      - compteurs du checkpoint == fichier
    #      - UNE SEULE empreinte de code
    #      - ZERO identite contradictoire
    #    (`portee_substitution.py` n'est PAS un critere de recette : il
    #     designe les lots ou sigma permute, ce qui reste vrai sur un
    #     fichier juste. Ce qu'on veut, c'est l'homogeneite.)
    # -----------------------------------------------------------------
    Etape 3 "Recette : compteurs, version unique, zero contradiction"
    python -u retirer_lots.py $dossier --verifier
    # Les ancres, cette fois contre le balayage. Sur un lot partiel, les
    # identites absentes sont annoncees comme telles et ne comptent pas
    # comme un accord -- un controle vide n'est pas un controle.
    python -u ancres_port.py $dossier
    Write-Host ""
    Write-Host ("  Session de {0:hh\:mm\:ss}." -f ((Get-Date) - $t0)) -ForegroundColor Green
    Write-Host "  Continuer :  .\run_propre.ps1 -SansTests" -ForegroundColor Green
    Write-Host "  Comparer les deux balayages, dans les DEUX sens, a la fin :" -ForegroundColor Green
    Write-Host "     python -u comparer_scans.py scan_wilson5 $dossier --sortie comparaison_w5_w6.json" -ForegroundColor Green
    Write-Host "  Quand le balayage annonce 0 lot a traiter :" -ForegroundColor Green
    Write-Host "               .\run_propre.ps1 -SansTests -ControleFinal" -ForegroundColor Green
    Write-Host ""
}
