# Microsoft Developer Studio Project File - Name="impex" - Package Owner=<4>
# Microsoft Developer Studio Generated Build File, Format Version 6.00
# ** NICHT BEARBEITEN **

# TARGTYPE "Win32 (x86) Static Library" 0x0104

CFG=impex - Win32 Debug
!MESSAGE Dies ist kein gültiges Makefile. Zum Erstellen dieses Projekts mit NMAKE
!MESSAGE verwenden Sie den Befehl "Makefile exportieren" und führen Sie den Befehl
!MESSAGE 
!MESSAGE NMAKE /f "impex.mak".
!MESSAGE 
!MESSAGE Sie können beim Ausführen von NMAKE eine Konfiguration angeben
!MESSAGE durch Definieren des Makros CFG in der Befehlszeile. Zum Beispiel:
!MESSAGE 
!MESSAGE NMAKE /f "impex.mak" CFG="impex - Win32 Debug"
!MESSAGE 
!MESSAGE Für die Konfiguration stehen zur Auswahl:
!MESSAGE 
!MESSAGE "impex - Win32 Release" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE "impex - Win32 Debug" (basierend auf  "Win32 (x86) Static Library")
!MESSAGE 

# Begin Project
# PROP AllowPerConfigDependencies 0
# PROP Scc_ProjName ""
# PROP Scc_LocalPath ""
CPP=cl.exe
RSC=rc.exe

!IF  "$(CFG)" == "impex - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release6.0"
# PROP Intermediate_Dir "Release6.0"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MD /W3 /GX /O2 /I "..\..\include" /I "..\..\..\jpeg-6b" /I "..\..\..\tiff-v3.5.7\libtiff" /D "NDEBUG" /D "WIN32" /D "_WINDOWS" /D "HasJPEG" /D "HasTIFF" /YX /FD /c
# ADD BASE RSC /l 0x407
# ADD RSC /l 0x407
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Release6.0\vigraimpex.lib"

!ELSEIF  "$(CFG)" == "impex - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /Z7 /Od /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /MDd /W3 /GX /Z7 /Od /I "..\..\..\tiffinclude" /I "..\..\..\jpeginclude" /I "..\..\include" /I "..\..\..\jpeg-6b" /I "..\..\..\tiff-v3.5.7\libtiff" /D "_DEBUG" /D "WIN32" /D "_WINDOWS" /D "HasJPEG" /D "HasTIFF" /YX /FD /c
# ADD BASE RSC /l 0x407
# ADD RSC /l 0x407
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"Debug6.0\vigraimpex.lib"

!ENDIF 

# Begin Target

# Name "impex - Win32 Release"
# Name "impex - Win32 Debug"
# Begin Source File

SOURCE=.\bmp.cxx
# End Source File
# Begin Source File

SOURCE=.\byteorder.cxx
# End Source File
# Begin Source File

SOURCE=.\codecmanager.cxx
# End Source File
# Begin Source File

SOURCE=.\gif.cxx
# End Source File
# Begin Source File

SOURCE=.\imageinfo.cxx
# End Source File
# Begin Source File

SOURCE=.\jpeg.cxx
# End Source File
# Begin Source File

SOURCE=.\png.cxx

!IF  "$(CFG)" == "impex - Win32 Release"

# SUBTRACT CPP /I "..\..\..\tiff-v3.5.7\libtiff"

!ELSEIF  "$(CFG)" == "impex - Win32 Debug"

!ENDIF 

# End Source File
# Begin Source File

SOURCE=.\pnm.cxx
# End Source File
# Begin Source File

SOURCE=.\sun.cxx
# End Source File
# Begin Source File

SOURCE=.\tiff.cxx
# End Source File
# Begin Source File

SOURCE=.\viff.cxx
# End Source File
# Begin Source File

SOURCE=.\void_vector.cxx
# End Source File
# End Target
# End Project
