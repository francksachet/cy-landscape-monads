# =====================================================================
#  relance_wilson4.ps1 -- chaine complete apres la correction du §5.23
#
#  A lancer depuis la racine du depot, sous PowerShell.
#  Chaque etape s'arrete si la precedente a echoue.
#
#  Duree totale : environ 2 h 15 (1 h de scan + 1 h d'equivariance).
# =====================================================================

$ErrorActionPreference = 'Stop'
$t0 = Get-Date

function Etape($n, $titre) {
    Write-Host ""
    Write-Host ("=" * 70) -ForegroundColor Cyan
    Write-Host "  ETAPE $n -- $titre" -ForegroundColor Cyan
    Write-Host ("=" * 70) -ForegroundColor Cyan
}

# ---------------------------------------------------------------------
Etape 1 "Tests de non-regression (34 tests, ~2 min)"
# Aucun scan ne doit partir sur un socle non verifie : c'est la seule
# etape dont l'echec doit tout arreter.
# ---------------------------------------------------------------------
python tests_regression.py
if ($LASTEXITCODE -ne 0) {
    Write-Host "  Tests en echec -- NE PAS lancer le scan." -ForegroundColor Red
    exit 1
}

# ---------------------------------------------------------------------
Etape 2 "Scan Wilson, generateur enumere (194 CICYs, ~1 h a 7 coeurs)"
# --output scan_wilson4 : on n'ecrase NI scan_wilson2 (seule trace des
# trois candidats produite avant la correction) NI scan_wilson3.
# ---------------------------------------------------------------------
python -u -m cy_landscape.main_optimized cicylist.txt -j 7 `
    --output scan_wilson4 --reset `
    --wilson wilson_cicys.json --max-charge 5 --n-random 5000 |
    Tee-Object -FilePath scan_wilson4.log

# ---------------------------------------------------------------------
Etape 3 "Audit (dedup, anomalie, coherence chi) -> results_clean.jsonl"
# --cicylist est necessaire : sans lui le controle d'anomalie n'a pas la
# geometrie et se tait (§5.22).
# ---------------------------------------------------------------------
python audit_results.py scan_wilson4 --cicylist cicylist.txt |
    Tee-Object -FilePath scan_wilson4_audit.log

# ---------------------------------------------------------------------
Etape 4 "Triage -> results_ranked.jsonl"
# ---------------------------------------------------------------------
python triage_clean.py scan_wilson4 |
    Tee-Object -FilePath scan_wilson4_triage.log

# ---------------------------------------------------------------------
Etape 5 "Equivariance, niveau des charges (~10 min)"
# CET ORDRE EST OBLIGATOIRE : equivariance.py produit le champ
# groupes_utiles dont equivariance_f.py se sert a l'etape suivante.
# ---------------------------------------------------------------------
python -u equivariance.py cicyquotients.m cicylist.txt scan_wilson4 |
    Tee-Object -FilePath scan_wilson4_equiv.log

# ---------------------------------------------------------------------
Etape 6 "Equivariance, niveau des polynomes (~1 h, mono-coeur)"
# ---------------------------------------------------------------------
python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson4 |
    Tee-Object -FilePath scan_wilson4_equiv_f.log

# ---------------------------------------------------------------------
Write-Host ""
Write-Host ("=" * 70) -ForegroundColor Green
Write-Host ("  TERMINE en {0:hh\:mm\:ss}" -f ((Get-Date) - $t0)) -ForegroundColor Green
Write-Host ("=" * 70) -ForegroundColor Green
Write-Host ""
Write-Host "  Le controle qui decide de tout -- les trois candidats sont-ils la ?"
Write-Host ""
Write-Host '    python -c "import json;s={r[''cicy''] for r in map(json.loads,open(''scan_wilson4/results.jsonl'',encoding=''utf-8''))};print({n:(n in s) for n in (6890,6947,6715)})"'
Write-Host ""
Write-Host "  Puis, sur les survivants :"
Write-Host "    python -u equivariance_f.py cicyquotients.m cicylist.txt scan_wilson4 --cicy 6890"
Write-Host ""
