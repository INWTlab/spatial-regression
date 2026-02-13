### Erweiterung für: Multivariate Kerndichteschätzung bei anonymisierten Geokoordinaten
* Wichtigstes Referenzpaper: https://academic.oup.com/jrsssa/article/180/1/161/7068203.
* Die Erweiterung ist das Einbeziehen eines Regressionsmodells bei der Schätzung der Dichte der Koordinaten.
* Die Vergleichsmethode ist im Paket KernelHeaping als Funktion dbivr implementiert: https://cran.r-project.org/web/packages/Kernelheaping/index.html.

### Nutzungshinweise
* Das Skript evaluation.R aus dem inst/scripts Ordner enthält den Code für den Start von Simulationsstudien, die die Erweiterung 
mit der dbivr Methode vergleichen.
* Das Skript nutzt dabei folgende Funktionalitäten, die als Funktionen in separaten Skripten bereitgestellt sind und per source-Befehl aufgerufen werden:
  * Generierung von simulierten Geokoordinaten-Datensätzen.
  * Schätzung der Dichte der Koordinaten mit der Erweiterung und der dbivr Methode.
  * Parametrisierte Funktion für die Durchführung einer Simulationsstudie und Berechnung der Vergleichsmetriken.
