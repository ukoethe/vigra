#define MyAppName "Vigra 1.4"
#define MyAppVersion "1.4.1a"
#define MyAppPublisher "Ullrich Köthe, University of Hamburg"
#define MyAppURL "http://kogs.informatik.uni-hamburg.de/~koethe/vigra/"
#define MySourcePath "C:\Users\koethe\src\vigra1.4.1a"
#define MyBinPath "C:\Users\koethe\src\vigra\current"
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
Source: {#MyBinPath}\include\external\fftw3.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\jconfig.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\jmorecfg.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\jpeglib.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\png.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\pngconf.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\tiff.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\tiffconf.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\tiffio.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\tiffvers.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\zconf.h; DestDir: {app}\include\external
Source: {#MyBinPath}\include\external\zlib.h; DestDir: {app}\include\external
Source: {#MyBinPath}\lib\FFTW3.lib; DestDir: {app}\lib
Source: {#MyBinPath}\lib\libjpeg.lib; DestDir: {app}\lib
Source: {#MyBinPath}\lib\libpng.lib; DestDir: {app}\lib
Source: {#MyBinPath}\lib\libtiff.lib; DestDir: {app}\lib
Source: {#MyBinPath}\lib\zlib.lib; DestDir: {app}\lib
Source: {#MyBinPath}\src\impex\vigraimpex.dll.lib; DestDir: {app}\lib
Source: {#MyBinPath}\src\impex\vigraimpex.dll; DestDir: {app}\bin
Source: {#MyBinPath}\src\examples\*.exe; DestDir: {app}\bin

[UninstallDelete]
Type: filesandordirs; Name: {app}\src\impex\Release
Type: filesandordirs; Name: {app}\src\impex\Debug
Type: files; Name: {app}\src\impex\*.exp
Type: files; Name: {app}\src\impex\*.pdb
Type: files; Name: {app}\src\impex\*.dll
Type: files; Name: {app}\src\impex\*.lib
Type: dirifempty; Name: {app}\src\impex

Type: filesandordirs; Name: {app}\src\examples\Release
Type: filesandordirs; Name: {app}\src\examples\Debug
Type: files; Name: {app}\src\examples\*.exe
Type: files; Name: {app}\src\examples\*.pdb
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
	     'You won''t be able to compile any programs using {#MyAppName}.');

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
