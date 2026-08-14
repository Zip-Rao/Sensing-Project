@ECHO OFF
REM Build the bilingual sqc documentation site on Windows.
REM   make html  -> English at build\html, Chinese at build\html\zh
pushd %~dp0

set SPHINXBUILD=sphinx-build
set SPHINXOPTS=-W --keep-going
set SOURCEDIR=source
set BUILDDIR=build

if "%1" == "clean" (
	rmdir /s /q "%BUILDDIR%" 2>nul
	del /q "%SOURCEDIR%\en\Simulation_sqc.ipynb" "%SOURCEDIR%\zh\Simulation_sqc.ipynb" 2>nul
	goto end
)

REM Copy the docs notebook (single source under source\_shared) into each
REM language tree (gitignored copies).
copy /y "%SOURCEDIR%\_shared\Simulation_sqc.ipynb" "%SOURCEDIR%\en\Simulation_sqc.ipynb" >nul
copy /y "%SOURCEDIR%\_shared\Simulation_sqc.ipynb" "%SOURCEDIR%\zh\Simulation_sqc.ipynb" >nul

if "%1" == "en" goto builden
if "%1" == "zh" goto buildzh

:buildall
%SPHINXBUILD% %SPHINXOPTS% -b html "%SOURCEDIR%\en" "%BUILDDIR%\html"
if errorlevel 1 goto end
%SPHINXBUILD% %SPHINXOPTS% -b html "%SOURCEDIR%\zh" "%BUILDDIR%\html\zh"
goto end

:builden
%SPHINXBUILD% %SPHINXOPTS% -b html "%SOURCEDIR%\en" "%BUILDDIR%\html"
goto end

:buildzh
%SPHINXBUILD% %SPHINXOPTS% -b html "%SOURCEDIR%\zh" "%BUILDDIR%\html\zh"
goto end

:end
popd
