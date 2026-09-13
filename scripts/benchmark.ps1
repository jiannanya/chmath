param(
    [Parameter(Mandatory=$true)][string]$Executable,
    [ValidateRange(16,16777216)][int]$Items=65536,
    [ValidateRange(0,62)][int]$Processor=0,
    [ValidateRange(1,20)][int]$Repeats=3,
    [switch]$TestSuite
)
$ErrorActionPreference='Stop'
$binary=(Resolve-Path -LiteralPath $Executable).Path
$affinity=[long]1 -shl $Processor
$allowed=[Diagnostics.Process]::GetCurrentProcess().ProcessorAffinity.ToInt64()
if (($allowed -band $affinity) -eq 0) {throw "Processor $Processor is not available in the current process affinity mask"}
"Binary: $binary"
"SHA256: $((Get-FileHash -LiteralPath $binary -Algorithm SHA256).Hash)"
"Logical processor=$Processor; items=$Items; process runs=$Repeats; UTC=$([DateTime]::UtcNow.ToString('o'))"
for($run=1;$run -le $Repeats;++$run) {
    $info=[Diagnostics.ProcessStartInfo]::new()
    $info.FileName=$binary
    $info.Arguments=if($TestSuite){''}else{[string]$Items}
    $info.UseShellExecute=$false
    $info.CreateNoWindow=$true
    $info.RedirectStandardOutput=$true
    $info.RedirectStandardError=$true
    $process=[Diagnostics.Process]::new()
    $process.StartInfo=$info
    try {
        if (!$process.Start()) {throw 'Benchmark process did not start'}
        $process.ProcessorAffinity=[IntPtr]$affinity
        $stdout=$process.StandardOutput.ReadToEnd()
        $stderr=$process.StandardError.ReadToEnd()
        $process.WaitForExit()
        if($process.ExitCode -ne 0){throw "Benchmark exited $($process.ExitCode): $stderr"}
        "Run $run"
        $stdout.TrimEnd()
        if($stderr){$stderr.TrimEnd()}
    } finally {$process.Dispose()}
}
