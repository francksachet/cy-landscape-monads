# =====================================================================
#  run_propre.ps1 -- le balayage d'equivariance recalcule A NEUF, dans
#  un seul etat du code.
#
#  POURQUOI A NEUF PLUTOT QU'UNE TROISIEME CHIRURGIE
#  --------------------------------------------------
#  `scan_wilson4` melange au moins trois versions du code. Deux
#  contaminations y ont ete trouvees le meme jour : le sigma du 5.34, et
#  des lignes « hors domaine » heritees de l'epoque ou `rank_c_max`
#  valait 1 -- 41 identites y portent A LA FOIS « hors domaine » et
#  « ok ». La seconde n'etait detectable que parce que `domaine_valide`
#  est une fonction pure qu'on peut rejouer pour rien. Les verdicts
#  couteux -- h0 equivariant, Hoppe, surjectivite -- issus des memes
#  vieux lots ne se revoient pas a ce prix, et rien ne dit qu'ils soient
#  epargnes.
#
#  RIEN N'EST DETRUIT. On ecrit dans scan_wilson5 ; scan_wilson4 reste
#  en place comme base de comparaison. Pas de --reset, donc pas d'etape
#  destructive du tout.
#
#  Duree : ~29 h 30 a 7 coeurs. Mesuree, pas estimee : la session du
#  17-18 aout a fait 21 968 realisations en 11 h 33, soit 1,89 s par
#  realisation, et il y en a 56 134 en tout. Fractionnable a volonte :
#  le checkpoint est ecrit apres chaque lot.
#
#      .\run_propre.ps1              # premiere fois
#      .\run_propre.ps1 -SansTests   # sessions suivantes
#      .\run_propre.ps1 -ControleFinal   # une fois TOUS les lots 'T' faits
# =====================================================================
param(
    [switch]$SansTests,
    [switch]$ControleFinal,
    [int]$Coeurs = 7
)

$ErrorActionPreference = 'Stop'
$dossier = 'scan_wilson5'
$source  = 'scan_wilson4'
$t0 = Get-Date

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
    Etape 1 "Tests de non-regression (45 tests, ~5 min)"
    python tests_regression.py
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  Tests en echec -- NE RIEN LANCER." -ForegroundColor Red
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
        --replier-orbites -j $Coeurs |
        Tee-Object -Append -FilePath scan_wilson5_equiv_f.log
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
    Write-Host ""
    Write-Host ("  Session de {0:hh\:mm\:ss}." -f ((Get-Date) - $t0)) -ForegroundColor Green
    Write-Host "  Continuer :  .\run_propre.ps1 -SansTests" -ForegroundColor Green
    Write-Host "  Quand le balayage annonce 0 lot a traiter :" -ForegroundColor Green
    Write-Host "               .\run_propre.ps1 -SansTests -ControleFinal" -ForegroundColor Green
    Write-Host ""
}
