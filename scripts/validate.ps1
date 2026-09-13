param([string]$Compiler = 'cl', [int]$Jobs = 8)
$ErrorActionPreference = 'Stop'
$projectRoot = Split-Path -Parent $PSScriptRoot
function Invoke-Checked {
    param([string]$Program, [string[]]$Arguments)
    & $Program @Arguments
    if ($LASTEXITCODE -ne 0) { throw "$Program failed with exit code $LASTEXITCODE" }
}
foreach ($config in @('Debug', 'Release', 'Scalar')) {
    $buildType = if ($config -eq 'Scalar') { 'Release' } else { $config }
    $simd = if ($config -eq 'Scalar') { 'OFF' } else { 'ON' }
    $compilerName = [IO.Path]::GetFileNameWithoutExtension($Compiler)
    $buildDir = Join-Path $projectRoot "build/validate-$compilerName-$config"
    Invoke-Checked 'cmake' @('-S', $projectRoot, '-B', $buildDir, '-G', 'Ninja', "-DCMAKE_CXX_COMPILER=$Compiler", "-DCMAKE_BUILD_TYPE=$buildType", "-DCHMATH_ENABLE_SIMD=$simd", '-DCHMATH_BUILD_TESTS=ON', '-DCHMATH_BUILD_EXAMPLES=ON', '-DCHMATH_BUILD_BENCHMARKS=ON')
    Invoke-Checked 'cmake' @('--build', $buildDir, '--parallel', "$Jobs")
    Invoke-Checked 'ctest' @('--test-dir', $buildDir, '--output-on-failure', '--parallel', "$Jobs")
}
