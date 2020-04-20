# Mallas contenidas en este directorio

## 1. SIANE_spain.geo
Malla generada a partir de los datos del Centro de Descargas del Ministerio de Transporte, http://centrodedescargas.cnig.es/CentroDescargas/
- En concreto, se usó el fichero "shapefile" `SIANE_CARTO_BASE_S_10M.zip`
- Se post-procesó con *QGis* y se exportó a `.geo` (ver notas abajo)
- Contiene fronteras entre comunidades autónomas
- De hecho, está separado entre CCAA, por lo que es difícil mallar todas españa

## 2. natural_spain.geo
Malla generada a partir de los datos de *Natural Earth*: https://www.naturalearthdata.com/downloads/110m-cultural-vectors/
- Se usó el fichero "shapefile" ne_10m_admin_0_map_subunits.zip
- Se filtró para extraer solamente a España
```
gr2ogr \
                -where "ADM0_A3 IN ('ESP')" \
                subunits.shp \
                ne_10m_admin_0_map_subunits.shp
```
- Se eliminaron manualmente islas, ciudades autónomas, isla de Peregil etc.
  - Activar edición (botón lápiz)
  - Menú Editar / Borrar parte
- Se post-procesó con *QGis* y se exportó a `.geo` (ver notas abajo)
- El plugin utiliza curvas Spline. Se convirtieron a rectas usando el
  script `geo_no_splines.py`, programado expresamente para ello:
    $ geo_no_splines fichero_entrada.geo > fichero_salida.geo

## Notas sobre proceso con Qgis
- Para exportar a .geo (*GMSH*), se utiliza el plugin
  *GMSH* (descargable directamente desde *QGis*)
- Éste necesita utilizar capas líneas
- Para usarlo, se generan líneas de frontera con
  "Vectorial/Herramientas de Geometría/Polígonos a Líneas"
