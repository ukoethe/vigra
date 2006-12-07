#define MyAppName "Vigra 1.5"
#define MyAppVersion "1.5.0"
#define MyAppPublisher "Ullrich Köthe, University of Hamburg"
#define MyAppURL "http://kogs.informatik.uni-hamburg.de/~koethe/vigra/"
#define MySourcePath "C:\Users\koethe\src\vigra1.5.0"
#define MyExternalPath "C:\Users\koethe\src\vigra\current"
#define RequiredVCVersion "7.1"

[Setup]
AppCopyright=Ullrich Köthe
AppName={#MyAppName}
AppVerName=Vigra {#MyAppVersion}
LicenseFile={#MySourcePath}\LICENSE.txt
ShowLanguageDialog=yes
AppSupportURL={#MyAppUrl}
AppVersion={#MyAppVersion}
DefaultDirName={reg:HKCU\Software\Microsoft\VisualStudio\7.1,VisualStudioProjectsLocation|{pf}}\{#MyAppName}
OutputBaseFilename=Setup-Vigra-{#MyAppVersion}
Compression=lzma
SolidCompression=true
AppID={#MyAppName}
OutputDir=C:\Users\koethe\src\vigra\current\WinSetup
InfoBeforeFile=

[Files]
Source: {#MySourcePath}\LICENSE.txt; DestDir: {app}
Source: {#MySourcePath}\README.txt; DestDir: {app}
Source: {#MySourcePath}\include\vigra\*.hxx; DestDir: {app}\include\vigra
Source: {#MySourcePath}\include\vigra\*.h; DestDir: {app}\include\vigra
Source: {#MySourcePath}\src\*.cxx; DestDir: {app}\src; Flags: recursesubdirs
Source: {#MySourcePath}\src\*.hxx; DestDir: {app}\src; Flags: recursesubdirs
Source: {#MySourcePath}\src\*.c; DestDir: {app}\src; Flags: recursesubdirs
Source: {#MySourcePath}\src\*.h; DestDir: {app}\src; Flags: recursesubdirs
Source: {#MySourcePath}\src\*.vcproj; DestDir: {app}\src; Flags: recursesubdirs
Source: {#MySourcePath}\src\*.sln; DestDir: {app}\src
Source: {#MySourcePath}\src\examples\*.gif; DestDir: {app}\src\examples
Source: {#MySourcePath}\test\*.cxx; DestDir: {app}\test; Flags: recursesubdirs
Source: {#MySourcePath}\test\*.hxx; DestDir: {app}\test; Flags: recursesubdirs
Source: {#MySourcePath}\test\*.gif; DestDir: {app}\test; Flags: recursesubdirs
Source: {#MySourcePath}\test\*.xv; DestDir: {app}\test; Flags: recursesubdirs
Source: {#MySourcePath}\test\*.vcproj; DestDir: {app}\test; Flags: recursesubdirs
Source: {#MySourcePath}\test\*.sln; DestDir: {app}\test
Source: {#MySourcePath}\test\testOrDelete.bat; DestDir: {app}\test
Source: {#MySourcePath}\doc\vigra\*; Excludes: CVS; DestDir: {app}\doc\vigra; Flags: recursesubdirs
Source: {#MySourcePath}\include\external\fftw3.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\jconfig.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\jmorecfg.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\jpeglib.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\png.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\pngconf.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\tiff.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\tiffconf.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\tiffio.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\tiffvers.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\zconf.h; DestDir: {app}\include\external
Source: {#MySourcePath}\include\external\zlib.h; DestDir: {app}\include\external
Source: {#MySourcePath}\lib\FFTW3.lib; DestDir: {app}\lib
Source: {#MySourcePath}\lib\libjpeg.lib; DestDir: {app}\lib
Source: {#MySourcePath}\lib\libpng.lib; DestDir: {app}\lib
Source: {#MySourcePath}\lib\libtiff.lib; DestDir: {app}\lib
Source: {#MySourcePath}\lib\zlib.lib; DestDir: {app}\lib
Source: {#MySourcePath}\src\impex\vigraimpex.dll.lib; DestDir: {app}\lib
Source: {#MySourcePath}\src\impex\vigraimpex.dll; DestDir: {app}\bin
Source: {#MySourcePath}\src\examples\*.exe; DestDir: {app}\bin

[UninstallDelete]
Type: filesandordirs; Name: {app}\src\impex\Release
Type: filesandordirs; Name: {app}\src\impex\Debug
Type: files; Name: {app}\src\impex\*.exp
Type: files; Name: {app}\src\impex\*.pdb
Type: files; Name: {app}\src\impex\*.dll
Type: files; Name: {app}\src\impex\*.lib
Type: files; Name: {app}\src\impex\*.ilk
Type: dirifempty; Name: {app}\src\impex

Type: filesandordirs; Name: {app}\src\examples\Release
Type: filesandordirs; Name: {app}\src\examples\Debug
Type: files; Name: {app}\src\examples\*.exe
Type: files; Name: {app}\src\examples\*.pdb
Type: files; Name: {app}\src\examples\*.ilk
Type: dirifempty; Name: {app}\src\examples

Type: files; Name: {app}\src\*.ncb
Type: files; Name: {app}\src\*.suo
Type: dirifempty; Name: {app}\src

Type: files; Name: {app}\lib\*.lib
Type: dirifempty; Name: {app}\lib

Type: filesandordirs; Name: {app}\test
Type: dirifempty; Name: {app}

[Code]
var
	VCPath: String;
	VCUserPath: String;

procedure CurStepChanged(CurStep: TSetupStep);
var
	Path: String;
begin
	if CurStep = ssPostInstall then
	begin
		RegQueryStringValue(HKCU, 'Environment', 'PATH', Path);
		Path := Path + ';' + ExpandConstant('{app}\bin');
		RegWriteStringValue(HKCU, 'Environment', 'PATH', Path);
	end;
end;

procedure CurUninstallStepChanged(CurUninstallStep: TUninstallStep);
var
	Path: String;
begin
	if CurUninstallStep = usPostUninstall then
	begin
		RegQueryStringValue(HKCU, 'Environment', 'PATH', Path);
		StringChange(Path,  ';' + ExpandConstant('{app}\bin'), '');
		RegWriteStringValue(HKCU, 'Environment', 'PATH', Path);
	end;
end;

procedure InitializeWizard();
begin
	if VCPath = '' then
		CreateOutputMsgPage(wpWelcome, 'Visual C++ {#RequiredVCVersion} not found.', '',
	     'Visual C++ {#RequiredVCVersion} could not be found on your computer. ' +
	     'You may need to re-compile {#MyAppName}, because the provided binaries ' +
	     'may be incompatible with your compiler.');

	if VCUserPath <> '' then
		WizardForm.DirEdit.Text := VCUserPath + '\{#MyAppName}'
	else
		WizardForm.DirEdit.Text := ExpandConstant('{pf}\{#MyAppName}');
end;

function InitializeSetup(): Boolean;
begin
	// Look for Visual Studio 7.1
	VCPath := '';
	VCUserPath := '';
	if RegQueryStringValue(HKLM, 'Software\Microsoft\VisualStudio\{#RequiredVCVersion}\Setup\VC', 'ProductDir', VCPath) then
		RegQueryStringValue(HKCU, 'Software\Microsoft\VisualStudio\{#RequiredVCVersion}','VisualStudioProjectsLocation', VCUserPath);
	Result := True;
end;
