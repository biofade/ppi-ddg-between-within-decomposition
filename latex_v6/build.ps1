# =====================================================================
#  build.ps1 -- compile the LaTeX conversion of v6
#  XeLaTeX + fontspec | equation env | vector PDF figures |
#  booktabs three-line tables | BibTeX (superscript numeric)
#
#  Pipe: xelatex -> bibtex -> reorder .bbl to v6 numbering -> xelatex x2
#  Run:  powershell -File build.ps1
# =====================================================================
$ErrorActionPreference = 'Continue'
$env:PATH = "F:\MiKTeX\miktex\bin\x64;$env:PATH"
$env:TMP  = "F:\dsh\cache\tmp"
$env:TEMP = "F:\dsh\cache\tmp"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location -LiteralPath $PSScriptRoot

Write-Host "=== manuscript_v6 (xelatex) ===" -ForegroundColor Cyan
Remove-Item *.aux,*.log,*.out,*.bbl,*.blg,*.pdf -EA 0
& xelatex -interaction=nonstopmode manuscript_v6.tex | Out-Null

Write-Host "=== bibtex ===" -ForegroundColor Cyan
& bibtex manuscript_v6 | Out-Null

Write-Host "=== reorder .bbl to v6 numbering ===" -ForegroundColor Cyan
python (Join-Path $repo 'scripts\v8_1\56_reorder_bbl.py') | Select-Object -First 1

Write-Host "=== manuscript_v6 (xelatex x2) ===" -ForegroundColor Cyan
& xelatex -interaction=nonstopmode manuscript_v6.tex | Out-Null
& xelatex -interaction=nonstopmode manuscript_v6.tex | Out-Null

Write-Host "=== supplementary_v6 (xelatex x3) ===" -ForegroundColor Cyan
1..3 | ForEach-Object { & xelatex -interaction=nonstopmode supplementary_v6.tex | Out-Null }

foreach ($f in 'manuscript_v6','supplementary_v6') {
    $log = Get-Content "$f.log" -Raw
    $pg  = (Select-String -Path "$f.log" -Pattern 'Output written on .* \((\d+) page' |
            Select-Object -Last 1)
    Write-Host ("  {0,-18} fatal={1}  citation_undef={2}  reference_undef={3}  pages={4}" -f `
        $f, ([regex]::Matches($log,'(?m)^!')).Count,
        ([regex]::Matches($log,'Citation .* undefined')).Count,
        ([regex]::Matches($log,'Reference .* undefined')).Count,
        $(if ($pg) { $pg.Matches.Groups[1].Value } else { '-' }))
}

Write-Host ""
Get-ChildItem *.pdf | ForEach-Object { "  {0,-28} {1,8:N0} B" -f $_.Name, $_.Length }
