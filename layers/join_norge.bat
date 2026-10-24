@echo off
rem Сборка архива модели сети Норвегии из частей
chcp 65001 >nul
copy /b elektro_topo_norge.zip.part00+elektro_topo_norge.zip.part01+elektro_topo_norge.zip.part02+elektro_topo_norge.zip.part03+elektro_topo_norge.zip.part04+elektro_topo_norge.zip.part05+elektro_topo_norge.zip.part06 elektro_topo_norge.zip
echo Готово: elektro_topo_norge.zip - распакуйте его, внутри elektro_topo_norge.gpkg
pause
