$ErrorActionPreference = "Stop"
$scriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$demoDir = $scriptDir

$PYTHON = "D:\conda\conda_envs\stimu\python.exe"
$TREE_SIM_ROOT = "C:\Users\intot\Desktop\molecular_clock\practice\seq_stimulation\reference"
$SEQ_SIM_ROOT  = "C:\Users\intot\Desktop\0811"

$env:PYTHONPATH = "$TREE_SIM_ROOT;$SEQ_SIM_ROOT;" + $env:PYTHONPATH

Write-Host "============================================"
Write-Host "  seq_sim_v2 Live Demo"
Write-Host "============================================"

Write-Host ""
Write-Host "[Step 1/4] Generating 5-tip tree (lognormal branch lengths) ..."
$treeFile = Join-Path $demoDir "demo_tree.nwk"
& $PYTHON -m tree_sim -n 5 -s 42 -c (Join-Path $demoDir "demo_lognormal.yaml") -o $treeFile

Write-Host "  Tree saved to: $treeFile"
Write-Host "  Tree content:"
Get-Content $treeFile
Write-Host ""

Write-Host "[Step 2/4] Simulating: equal base frequency (A=C=G=T=0.25) ..."
$outDir = Join-Path $demoDir "out_equal"
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
& $PYTHON -m seq_sim_v2 -m HKY -l 2000 -d 0.5 -z 42 -o f -y $outDir -q $treeFile
Write-Host "  Output: $outDir"

Write-Host "[Step 3/4] Simulating: high GC (A=0.1 C=0.4 G=0.4 T=0.1) ..."
$outDir = Join-Path $demoDir "out_gc"
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
& $PYTHON -m seq_sim_v2 -m HKY -l 2000 -d 0.5 -z 42 -f 0.1 0.4 0.4 0.1 -o f -y $outDir -q $treeFile
Write-Host "  Output: $outDir"

Write-Host "[Step 4/4] Simulating: high AT (A=0.4 C=0.1 G=0.1 T=0.4) ..."
$outDir = Join-Path $demoDir "out_at"
New-Item -ItemType Directory -Path $outDir -Force | Out-Null
& $PYTHON -m seq_sim_v2 -m HKY -l 2000 -d 0.5 -z 42 -f 0.4 0.1 0.1 0.4 -o f -y $outDir -q $treeFile
Write-Host "  Output: $outDir"

Write-Host ""
Write-Host "============================================"
Write-Host "  All simulations complete!"
Write-Host ""
Write-Host "  Equal freq: $demoDir\out_equal\"
Write-Host "  High GC:    $demoDir\out_gc\"
Write-Host "  High AT:    $demoDir\out_at\"
Write-Host ""
Write-Host "  Next: python count_freqs.py <fasta_file>"
Write-Host "============================================"
