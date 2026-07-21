# Set variables
$RscriptPath = "C:\Program Files\R\R-4.6.1\bin\x64\Rscript.exe"
$scriptPath = "long_sims Stan.R"
$logFolder = ".\logs"

# Create log folder if it doesn't exist
if (-not (Test-Path -Path $logFolder)) {
    New-Item -ItemType Directory -Path $logFolder
}

# Launch 10 batch jobs in parallel, each logging its output
for ($i = 1; $i -le 10; $i++) {
    $logFile = Join-Path $logFolder ("batch_$i.log")
    $errorFile = Join-Path $logFolder ("batch_${i}_error.log")

    Start-Process $RscriptPath `
        -ArgumentList "`"$scriptPath`" $i" `
        -RedirectStandardOutput $logFile `
        -RedirectStandardError $errorFile `
        -NoNewWindow
}

Write-Host "All batch jobs started!"
