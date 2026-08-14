$ErrorActionPreference = 'Stop'

Push-Location $PSScriptRoot
try {
    & latexmk -pdf main.tex
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }

    & latexmk -pdf supplemental.tex
    if ($LASTEXITCODE -ne 0) { exit $LASTEXITCODE }
}
finally {
    Pop-Location
}
