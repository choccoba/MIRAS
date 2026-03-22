# MIRAS Project Structure Setup
# Run from: G:\내 드라이브\MIRAS
# Usage: cd "G:\내 드라이브\MIRAS" && .\setup_MIRAS.ps1

$ROOT = Split-Path -Parent $MyInvocation.MyCommand.Path
Write-Host "MIRAS Setup: $ROOT" -ForegroundColor Cyan

$folders = @(
    "01_Step1_GEO_DEG\raw_data",
    "01_Step1_GEO_DEG\results",
    "01_Step1_GEO_DEG\figures",
    "02_Step2_miRNA_Mapping\results",
    "02_Step2_miRNA_Mapping\figures",
    "03_Step3_MetaAnalysis\literature_db",
    "03_Step3_MetaAnalysis\results",
    "03_Step3_MetaAnalysis\figures",
    "04_Step4_Integration\network",
    "04_Step4_Integration\results",
    "04_Step4_Integration\figures",
    "05_Reference_DBs\HMDD_v4",
    "05_Reference_DBs\miRBase",
    "06_v18_Connection\anchor_genes",
    "06_v18_Connection\overlap_analysis",
    "07_Disease_Map\drafts",
    "07_Disease_Map\figures",
    "08_Manuscripts\draft_v1",
    "99_Logs\session_notes",
    "99_Logs\claude_sessions"
)

foreach ($f in $folders) {
    $path = Join-Path $ROOT $f
    New-Item -ItemType Directory -Force -Path $path | Out-Null
    New-Item -ItemType File -Force -Path "$path\.gitkeep" | Out-Null
    Write-Host "  [+] $f" -ForegroundColor Gray
}

Write-Host ""
Write-Host "Committing folder structure to GitHub..." -ForegroundColor Cyan
Set-Location $ROOT
git add .
git commit -m "feat: initialize MIRAS analysis pipeline folder structure"
git push origin main

Write-Host ""
Write-Host "MIRAS setup complete!" -ForegroundColor Green
Write-Host "GitHub  : https://github.com/choccoba/MIRAS" -ForegroundColor White
Write-Host "Local   : $ROOT" -ForegroundColor White
