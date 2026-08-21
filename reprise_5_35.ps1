# =====================================================================
#  reprise_5_35.ps1 -- recalcul des lots fausses par le defaut du 5.34
#
#  A lancer autant de fois qu'on veut. Chaque passage reprend ou le
#  precedent s'est arrete. Ctrl-C est une facon NORMALE de s'arreter.
#
#  POURQUOI PASSER PAR CE SCRIPT PLUTOT QUE DE RETAPER LA COMMANDE
#  ---------------------------------------------------------------
#  L'empreinte du checkpoint est calculee sur (--cicy, --replier-orbites,
#  --controle-orbites, --taille-lot). Si l'une de ces quatre options
#  differe d'une session a l'autre -- un --replier-orbites oublie suffit --
#  `equivariance_f.py` declare le checkpoint inutilisable, EFFACE le JSONL
#  et repart de zero. Soit, ici, 160 Mo et plusieurs dizaines d'heures.
#
#  Ce script fige la ligne. Ne pas la modifier : la relancer telle quelle
#  est ce qui rend l'interruption sans consequence.
#
#  Duree totale estimee : ~11 h a 7 coeurs (21 968 realisations), a
#  prendre comme un ordre de grandeur -- les lots de #480 sont plus
#  lourds que la moyenne.
# =====================================================================
param(
    [switch]$SansTests,   # sauter les tests de non-regression
    [int]$Coeurs = 7
)

$ErrorActionPreference = 'Stop'
$dossier = 'scan_wilson4'
$t0 = Get-Date

function Etape($n, $titre) {
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "  ETAPE $n -- $titre" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
}

# ---------------------------------------------------------------------
# 1. Le socle. Aucun calcul ne part sur un socle non verifie.
#    -SansTests aux sessions suivantes : rien n'a change entre-temps.
# ---------------------------------------------------------------------
if (-not $SansTests) {
    Etape 1 "Tests de non-regression (42 tests, ~5 min)"
    python tests_regression.py
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  Tests en echec -- NE RIEN LANCER." -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "  (tests sautes : -SansTests)" -ForegroundColor DarkGray
}

# ---------------------------------------------------------------------
# 2. Le retrait. Ne s'execute qu'une fois : le script pose un marqueur et
#    refuse de se rejouer. C'est necessaire -- rejoue APRES la relance,
#    il reproposerait de retirer les lots qui viennent d'etre recalcules.
# ---------------------------------------------------------------------
$marqueur = Join-Path $dossier 'retrait_5_35_applique.json'
if (-not (Test-Path $marqueur)) {
    Etape 2 "Retrait des lots fausses (12 627 lignes, ~4 min)"
    Write-Host "  Cette etape n'est pas reprenable. Elle est courte." -ForegroundColor Yellow
    Write-Host "  Si elle est coupee pendant le remplacement des fichiers :" -ForegroundColor Yellow
    Write-Host "    restaurer les deux fichiers depuis ${dossier}_avant_5.35\" -ForegroundColor Yellow
    Write-Host "    puis relancer avec --refaire." -ForegroundColor Yellow
    python -u retirer_lots.py $dossier --appliquer
    if ($LASTEXITCODE -ne 0) {
        Write-Host "  Retrait en echec -- ne pas enchainer." -ForegroundColor Red
        exit 1
    }
} else {
    Write-Host "  (retrait deja applique)" -ForegroundColor DarkGray
}

# ---------------------------------------------------------------------
# 3. Le calcul. LA ligne, invariante d'une session a l'autre.
#    Pas de --max-realisations : sans plafond.
#    Le log est en -Append pour ne pas perdre les sessions precedentes.
#    Le checkpoint est ecrit apres CHAQUE lot ; une coupure coute au pire
#    le lot en cours.
# ---------------------------------------------------------------------
Etape 3 "Equivariance de f, sans plafond (Ctrl-C quand vous voulez)"
try {
    python -u equivariance_f.py cicyquotients.m cicylist.txt $dossier `
        --replier-orbites -j $Coeurs |
        Tee-Object -Append -FilePath scan_wilson4_equiv_f_5.35.log
} finally {
    # ---------------------------------------------------------------
    # 4. Controle d'integrite, execute MEME apres un Ctrl-C : les
    #    compteurs du checkpoint doivent decrire le fichier. Un
    #    desaccord signale un lot enregistre dont les lignes manquent.
    # ---------------------------------------------------------------
    Etape 4 "Controle d'integrite"
    python -u retirer_lots.py $dossier --verifier
    Write-Host ""
    Write-Host ("  Session de {0:hh\:mm\:ss}." -f ((Get-Date) - $t0)) -ForegroundColor Green
    Write-Host "  Relancer CE MEME script pour continuer :" -ForegroundColor Green
    Write-Host "    .\reprise_5_35.ps1 -SansTests" -ForegroundColor Green
    Write-Host ""
}
