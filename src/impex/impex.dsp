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
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
# ADD BASE CPP /nologo /W3 /GX /O2 /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /YX /FD /c
# ADD CPP /nologo /W3 /GX /O2 /I "..\..\include" /I "..\..\..\jpeginclude" /I "..\..\..\tiffinclude" /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "HasJPEG" /D "HasTIFF" /YX /FD /c
# ADD BASE RSC /l 0x407
# ADD RSC /l 0x407
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\lib\windows\vigraimpex.lib"

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
# ADD CPP /nologo /W3 /GX /Z7 /Od /I "..\..\..\tiffinclude" /I "..\..\..\jpeginclude" /I "..\..\..\jpeg-6b" /I "..\..\include" /D "WIN32" /D "_DEBUG" /D "_WINDOWS" /D "HasTIFF" /D "HasJPEG" /YX /FD /c
# ADD BASE RSC /l 0x407
# ADD RSC /l 0x407
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
LIB32=link.exe -lib
# ADD BASE LIB32 /nologo
# ADD LIB32 /nologo /out:"..\..\lib\windows\vigraimpex.lib"

!ENDIF 

# Begin Target

# Name "impex - Win32 Release"
# Name "impex - Win32 Debug"
# Begin Source File

SOURCE=.\block_read.c
# End Source File
# Begin Source File

SOURCE=.\bmp.c
# End Source File
# Begin Source File

SOURCE=.\colors.c
# End Source File
# Begin Source File

SOURCE=.\cr_image.c
# End Source File
# Begin Source File

SOURCE=.\error.c
# End Source File
# Begin Source File

SOURCE=.\freeimage.c
# End Source File
# Begin Source File

SOURCE=.\fullpath.c
# End Source File
# Begin Source File

SOURCE=.\gif.c
# End Source File
# Begin Source File

SOURCE=.\image.c
# End Source File
# Begin Source File

SOURCE=.\imagesize.c
# End Source File
# Begin Source File

SOURCE=.\impex.cxx
# End Source File
# Begin Source File

SOURCE=.\jpeg.c
# End Source File
# Begin Source File

SOURCE=.\machorder.c
# End Source File
# Begin Source File

SOURCE=.\machtype.c
# End Source File
# Begin Source File

SOURCE=.\nt.c
# End Source File
# Begin Source File

SOURCE=.\order.c
# End Source File
# Begin Source File

SOURCE=.\pnm.c
# End Source File
# Begin Source File

SOURCE=.\quantize.c
# End Source File
# Begin Source File

SOURCE=.\readheader.c
# End Source File
# Begin Source File

SOURCE=.\readimage.c
# End Source File
# Begin Source File

SOURCE=.\sun.c
# End Source File
# Begin Source File

SOURCE=.\tiff.cxx
# End Source File
# Begin Source File

SOURCE=.\utility.c
# End Source File
# Begin Source File

SOURCE=.\utils.c
# End Source File
# Begin Source File

SOURCE=.\viff.c
# End Source File
# Begin Source File

SOURCE=.\writeimage.c
# End Source File
# End Target
# End Project
