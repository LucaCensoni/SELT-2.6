
import sys
import os
import numpy as np
from qgis.core import (QgsApplication, QgsRasterLayer, QgsProject, QgsVectorLayer,
QgsWkbTypes, QgsProcessingException, QgsGraduatedSymbolRenderer,
QgsRendererRange, QgsSymbol, QgsGradientColorRamp,
QgsTextFormat, QgsTextBufferSettings, QgsPalLayerSettings, QgsVectorLayerSimpleLabeling)
from qgis.analysis import QgsRasterCalculator, QgsRasterCalculatorEntry
from osgeo import gdal
from qgis.PyQt.QtGui import QColor, QFont
from qgis.utils import iface

#CARICARE IL DTM MANUALMENTE

###     DATI    ###
#VASCA
volume_vasca = 30000.00
volume_min = 25000.00
volume_max = 35000.00
DHSu = 2 #passo con il quale si alza lo sbarramento

#PARAMETRI DI CALCOLO
passo_di_calcolo = 50
EPSG = 32632
valore = 1000000 #valore iniziale per la ricerca dei corsi d'acqua'
#RICERCA DELL'ALTEZZA DI SBARRAMENTO
errore = float('inf')
deltaerrore = 5000 #errore = volume calcolato - volume_vasca
ncicli = 30 #numero di cicli oltre il quale la ricerca si interrompe

#CARTELLA DOVE FAR AVVENIRE IL CALCOLO
Path_cartella_lavoro = "D:/Tesi Magistrale/TEST_1" #DOVE SALVARE I DATI DELL'ALGORITMO


#CREAZIONE DELLE CARTELLE DI LAVORO

directory_Area = f"{Path_cartella_lavoro}/Area"
if not os.path.exists(directory_Area):
    os.makedirs(directory_Area)
else:
    print(f"La cartella Area esiste già.")
directory_Bacino = f"{Path_cartella_lavoro}/Bacino"
if not os.path.exists(directory_Bacino):
    os.makedirs(directory_Bacino)
else:
    print(f"La cartella Bacino esiste già.")
directory_Vasca = f"{Path_cartella_lavoro}/Vasca"
if not os.path.exists(directory_Vasca):
    os.makedirs(directory_Vasca)
else:
    print(f"La cartella Vasca esiste già.")
directory_Vect = f"{Path_cartella_lavoro}/Vect"
if not os.path.exists(directory_Vect):
    os.makedirs(directory_Vect)
else:
    print(f"La cartella Vect esiste già.")
directory_Polig = f"{Path_cartella_lavoro}/Polig"
if not os.path.exists(directory_Polig):
    os.makedirs(directory_Polig)
else:
    print(f"La cartella Polig esiste già.")
directory_Extract = f"{Path_cartella_lavoro}/Extract"
if not os.path.exists(directory_Extract):
    os.makedirs(directory_Extract)
else:
    print(f"La cartella Extract esiste già.")
directory_Zona = f"{Path_cartella_lavoro}/Zona"
if not os.path.exists(directory_Zona):
    os.makedirs(directory_Zona)
else:
    print(f"La cartella Zona esiste già.")
directory_Diff = f"{Path_cartella_lavoro}/Diff"
if not os.path.exists(directory_Diff):
    os.makedirs(directory_Diff)
else:
    print(f"La cartella Diff esiste già.")

#DTM
# DTM_file_path = f"{Path_cartella_lavoro}/DTM/{Nome_file_DTM}"
# DTM = QgsRasterLayer(DTM_file_path,"DTM")

# if not DTM.isValid():
#     print("Layer DTM non valido")
# else:
#     QgsProject.instance().addMapLayer(DTM)

#Algoritmo Fill & Sink

processing.run("sagang:fillsinksplanchondarboux2001", {'DEM':'DTM',
'RESULT':f"{Path_cartella_lavoro}/DTM_FILL.sdat",
'MINSLOPE':0.01})
DTM_FILL_file_path = f"{Path_cartella_lavoro}/DTM_FILL.sdat"
DTM_FILL = QgsRasterLayer(DTM_FILL_file_path, "DTM_FILL")

if not DTM_FILL.isValid():
    print("Layer DTM_FILL non valido")
else:
    QgsProject.instance().addMapLayer(DTM_FILL)


print("############################")
print("CARICAMENTO LAYER DTM e FILL & SINK EFFETTUATO!!")
print("############################")

#Algoritmo Catchmentarea

processing.run("sagang:catchmentarea", {'ELEVATION':'DTM_FILL',
'SINKROUTE':None,
'WEIGHTS':None,
'FLOW':f"{Path_cartella_lavoro}/Area_inflow_D8.sdat",
'VAL_INPUT':None,
'VAL_MEAN':'TEMPORARY_OUTPUT',
'ACCU_MATERIAL':None,
'ACCU_TARGET':'DTM_FILL',
'ACCU_TOTAL':'TEMPORARY_OUTPUT',
'ACCU_LEFT':'TEMPORARY_OUTPUT',
'ACCU_RIGHT':'TEMPORARY_OUTPUT',
'STEP':1,
'FLOW_UNIT':1,
'FLOW_LENGTH':f"{Path_cartella_lavoro}/path_length.sdat",
'LINEAR_VAL':None,
'LINEAR_DIR':None,
'WEIGHT_LOSS':'TEMPORARY_OUTPUT',
'METHOD':0,
'LINEAR_DO':False,
'LINEAR_MIN':500,
'CONVERGENCE':1.1,
'MFD_CONTOUR':False,
'NO_NEGATIVES':True})

Area_inflow_D8_file_path = f"{Path_cartella_lavoro}/Area_inflow_D8.sdat"
Area_inflow_D8 = QgsRasterLayer(Area_inflow_D8_file_path, "Flow_D8")
if not Area_inflow_D8.isValid():
    print("Layer Area_infow_D8 non valido")
else:
    QgsProject.instance().addMapLayer(Area_inflow_D8)
    
Path_length_file_path = f"{Path_cartella_lavoro}/path_length.sdat"
Path_length = QgsRasterLayer(Path_length_file_path, "Path_length")
if not Path_length.isValid():
    print("Layer Path_length non valido")
else:
    QgsProject.instance().addMapLayer(Path_length)


print("############################")
print("Catchmentarea completato con caricamento Flow_D8 e Path Length!!")
print("############################")

#ALGORITMO CHANNEL NETWORK

processing.run("sagang:channelnetwork",
{'ELEVATION':'DTM_FILL',
'SINKROUTE':None,
'CHNLNTWRK':f"{Path_cartella_lavoro}/ch_net_{valore}.sdat",
'CHNLROUTE':f"{Path_cartella_lavoro}/ch_dir_{valore}.sdat",
'SHAPES':f"{Path_cartella_lavoro}/ch_net.shp",
'INIT_GRID':f"{Path_cartella_lavoro}/Area_inflow_D8.sdat",
'INIT_METHOD':2,
'INIT_VALUE': valore,
'DIV_GRID':None,
'DIV_CELLS':1,
'TRACE_WEIGHT':None,
'MINLEN':0})

Ch_net_file_path = f'C:/Users/Luca.Censoni/Desktop/Modifica raster/ch_net_{valore}.sdat'
Ch_net = QgsRasterLayer(Ch_net_file_path, f"Ch_net_{valore}")
if not Ch_net.isValid():
    print(f"Layer Ch_net_{valore} non valido")
else:
    QgsProject.instance().addMapLayer(Ch_net)

Ch_dir_file_path = f'C:/Users/Luca.Censoni/Desktop/Modifica raster/ch_dir_{valore}.sdat'
Ch_dir = QgsRasterLayer(Ch_dir_file_path, f"Ch_dir_{valore}")
if not Ch_dir.isValid():
    print(f"Layer Ch_dir_{valore} non valido")
else:
    QgsProject.instance().addMapLayer(Ch_dir)

Ch_net_shp_file_path = f"{Path_cartella_lavoro}/ch_net.shp"
Ch_net_shp = QgsVectorLayer(Ch_net_shp_file_path, "Ch_net_shp")
(Ch_net_shp_file_path, "Ch_net_shp")
if not Ch_net_shp.isValid():
    print("Layer Ch_net_shp non valido")
else:
    QgsProject.instance().addMapLayer(Ch_net_shp)

# Stile del layer dei fiumi
stile_linea = QgsLineSymbol.createSimple({'color': 'blu', 'width': '0.3'})
render_linea = QgsSingleSymbolRenderer(stile_linea)
Ch_net_shp.setRenderer(render_linea)
QgsProject.instance().addMapLayer(Ch_net_shp)

print("############################")
print("CHANNEL NETWORK COMPLETATO!!")
print("############################")
#SEZIONI DI CHIUSURA

# Creazione di un nuovo layer per i punti
layer_sez_chiusura = QgsVectorLayer(f"Point?crs=EPSG:{EPSG}", "Sez_chiusura", "memory")
layer_sez_chiusura_provider = layer_sez_chiusura.dataProvider()

# Creazione dei punti lungo le geometrie del layer dei fiumi
distanza = QgsDistanceArea()
distanza.setSourceCrs(Ch_net_shp.crs(), QgsProject.instance().transformContext())
distanza.setEllipsoid(Ch_net_shp.crs().ellipsoidAcronym())

for feature in Ch_net_shp.getFeatures():
    geom = feature.geometry()
    length = distanza.measureLength(geom)

    current_distance = 0
    while current_distance < length:
        point = geom.interpolate(current_distance)
        point_feature = QgsFeature()
        point_feature.setGeometry(point)
        layer_sez_chiusura_provider.addFeature(point_feature)
        current_distance += passo_di_calcolo

# Aggiunta dell'elevazione dal DTM ai punti
layer_sez_chiusura.startEditing()
layer_sez_chiusura.dataProvider().addAttributes([
    QgsField('x_coord', QVariant.Double),
    QgsField('y_coord', QVariant.Double),
    QgsField('elevation', QVariant.Double),
    QgsField('zs',QVariant.Double),
    QgsField('Area Bac',QVariant.Double),
    QgsField('Volume Bac', QVariant.Double),
    QgsField('HS',QVariant.Double)
])
layer_sez_chiusura.updateFields()

#Ciclo per aggiungere le cooordinate e quota

for feature in layer_sez_chiusura.getFeatures():
    geom = feature.geometry()
    if geom.type() == QgsWkbTypes.PointGeometry:
        x = geom.asPoint().x()
        y = geom.asPoint().y()
        feature['x_coord'] = x
        feature['y_coord'] = y

        # Identifica l'elevazione dal DTM per il punto
        intersected_feature = DTM_FILL.dataProvider().identify(geom.asPoint(), QgsRaster.IdentifyFormatValue).results()
        if intersected_feature:
            value = intersected_feature[1]
            if value is not None:
                feature['elevation'] = round(value, 3)

        layer_sez_chiusura.updateFeature(feature)

layer_sez_chiusura.commitChanges()

# Salvataggio del layer vettoriale dei punti
output_points_path = f"{Path_cartella_lavoro}/Punti.shp"
QgsVectorFileWriter.writeAsVectorFormat(layer_sez_chiusura, output_points_path, "UTF-8", layer_sez_chiusura.crs(), "ESRI Shapefile")

# Aggiunta del layer vettoriale dei punti al progetto QGIS
QgsProject.instance().addMapLayer(QgsVectorLayer(output_points_path, "Sez_Chiusura", "ogr"))

count = layer_sez_chiusura.featureCount()
print(f"totale punti = {count}")


#ALGORITMO PRE OGNI PUNTO INDIVIDUATO LUNGO I FIUMI


def apply_upslope_area_to_points(layer_sez_chiusura, DHSu):
        counter = 1
        
        if layer_sez_chiusura.geometryType() == QgsWkbTypes.PointGeometry:
            for feature in layer_sez_chiusura.getFeatures():
                x_ = feature.geometry()
                attrs = feature.attributes()
                x = attrs[0]
                y = attrs[1]
                Q = attrs[2]
                print(f"{(counter/count) *100}%")
                print(f"Feature ID {feature.id()} - X: {x}, Y: {y}, Q: {Q}")
            
               # Crea un feedback per mostrare il progresso
                feedback = QgsProcessingFeedback()
                
                #UPSLOPE AREA
                processing.run("sagang:upslopearea", {'TARGET':None,
                'TARGET_PT_X':x,
                'TARGET_PT_Y':y,
                'ELEVATION':'DTM_FILL',
                'SINKROUTE':None,
                'AREA':f"{Path_cartella_lavoro}/Area/Area_{counter}.sdat",
                'METHOD':2,
                'CONVERGE':1.1,
                'MFD_CONTOUR':False},
                feedback=feedback)
    
                # Ottieni il percorso del risultato
                area_file_path = f"{Path_cartella_lavoro}/Area/Area_{counter}.sdat"
                area_layer = QgsRasterLayer(area_file_path, f'Area_{counter}')

                if not area_layer.isValid():
                    print(f"Layer Area_{counter} non valido")
                else:
                    QgsProject.instance().addMapLayer(area_layer)
                    
                #POLYGONIZE (RASTER TO VECTOR)
                processing.run("gdal:polygonize", {'INPUT':f"{Path_cartella_lavoro}/Area/Area_{counter}.sdat",
                'BAND':1,
                'FIELD':'DN',
                'EIGHT_CONNECTEDNESS':False,
                'EXTRA':'',
                'OUTPUT':f"{Path_cartella_lavoro}/Vect/Vect_{counter}.shp"})
                Vect_file_path = f"{Path_cartella_lavoro}/Vect/Vect_{counter}.shp"
                QgsProject.instance().addMapLayer(QgsVectorLayer(Vect_file_path, f"Vect_{counter}", "ogr"))
                
                #CLIP I DTM_FILL
                processing.run("sagang:cliprasterwithpolygon", {'INPUT':[f'{Path_cartella_lavoro}/DTM_FILL.sdat'],
                'OUTPUT':f'{Path_cartella_lavoro}/Bacino/Bacino_{counter}.sdat',
                'POLYGONS':f'{Path_cartella_lavoro}/Vect/Vect_{counter}.shp',
                'EXTENT':1})
                
                # carica il risultato
                bacino_file_path = f"{Path_cartella_lavoro}/Bacino/Bacino_{counter}.sdat"
                bacino_layer = QgsRasterLayer(bacino_file_path, f'Bacino_{counter}')

                if not bacino_layer.isValid():
                    print(f"Layer Bacino_{counter} non valido")
                else:
                    QgsProject.instance().addMapLayer(bacino_layer)
                
                #creazione delle liste
                volumi_invaso = []
                layer_da_rimuovere = []
                layers_1 = []
                errore = float('inf')
                #print(DHSu)
                DHS = DHSu
                #print(f'DHS = {DHS}')
                deltah = DHSu
                #print(deltah)
                counter2 = 1
                counter3 = 1
                while errore > deltaerrore and counter2 < ncicli:
                #CALCOLO DEL PERIMETRO 
                    zs = Q +deltah #quota sommitale sbarramento
                    #print(f"zs={zs}")
                    espressione = f'(Bacino_{counter}@1) < {zs}'
                
                    raster_lavorato = QgsProject.instance().mapLayersByName(f'Bacino_{counter}')[0]
                    #DEFINISCI VOCI DEL RASTER DA CALCOLARE 
                    raster_lavorato_entry = QgsRasterCalculatorEntry()
                    raster_lavorato_entry.ref = f'Bacino_{counter}@1'
                    raster_lavorato_entry.raster = raster_lavorato
                    raster_lavorato_entry.bandNumber = 1
                
                    #QgsRASTERCALCULATOR CON ESPRESSIONE VARIABILE
                    calc = QgsRasterCalculator(espressione,f'{Path_cartella_lavoro}/Vasca/Vasca_{counter}_{zs}.tif',
                        'GTiff', raster_lavorato.extent(), raster_lavorato.width(),
                        raster_lavorato.height(), [raster_lavorato_entry])
                    
                #ESEGUI IL CALCOLO
                    calc.processCalculation()
                
                    output_raster_lavorato_path = f'{Path_cartella_lavoro}/Vasca/Vasca_{counter}_{zs}.tif'
                    output_raster_lavorato_name = f'Vasca_{counter}_{zs}'
                    output_raster_lavorato = QgsRasterLayer(output_raster_lavorato_path, output_raster_lavorato_name, 'gdal')
                    QgsProject.instance().addMapLayer(output_raster_lavorato)
            
            
                    processing.run("native:rastercalc", {'LAYERS':['DTM_FILL',
                    f'Vasca_{counter}_{zs}'],
                    'EXPRESSION':f'"Vasca_{counter}_{zs}@1" * "DTM_FILL@1"',
                    'EXTENT':None,
                    'CELL_SIZE':None,
                    'CRS':QgsCoordinateReferenceSystem(f'EPSG:{EPSG}'),
                    'OUTPUT':f'{Path_cartella_lavoro}/Vasca/Invaso_{counter}_{zs}.tif'})
            
                    invaso_file_path = f"{Path_cartella_lavoro}/Vasca/Invaso_{counter}_{zs}.tif"
                    invaso_layer = QgsRasterLayer(invaso_file_path, f'Invaso_{counter}_{zs}')
                    QgsProject.instance().addMapLayer(invaso_layer)
                
                    #CALCOLO DATI DELLA AREA
                    if not invaso_layer.isValid():
                        print(f"Invaso_{counter} è un file non è valido!")
                    else:
                        # Ottieni il data provider
                        provider = invaso_layer.dataProvider()
                    # Ottieni le dimensioni del raster
                    rows = provider.ySize()
                    cols = provider.xSize()
                
                    # Ottieni il blocco di dati raster
                    block = provider.block(1, invaso_layer.extent(), cols, rows)
                
                    # Conta le celle attive
                    active_cells = 0
                    no_data_value = provider.sourceNoDataValue(1)
                    for row in range(rows):
                        for col in range(cols):
                            value = block.value(row, col)
                            if value != no_data_value and value > 0:
                                active_cells += 1
                    # Ottieni la risoluzione spaziale del raster
                    cell_width = invaso_layer.rasterUnitsPerPixelX()
                    cell_height = invaso_layer.rasterUnitsPerPixelY()
                
                
                    layer_sez_chiusura.startEditing()
                    # print(f"Celle attive: {active_cells}")
                    area_cella = cell_width * cell_height
                    # print(f"Dim. cella = {area_cella}m^2")
                    area_totale = area_cella*active_cells
                    area_formattata = format(area_totale, ".1f")
                    # print(f"Area totale: {area_totale} m^2")
                
                    #CALCOLO DEL VOLUME
                    #Nome del layer DTM
                    nome_layer_dtm = f'Invaso_{counter}_{zs}'

                    # Recupera il layer raster del DTM
                    raster_dtm = QgsProject.instance().mapLayersByName(nome_layer_dtm)[0]

                    # Definisci l'espressione per il calcolo della differenza tra il DTM e la quota orizzontale
                    espressione_differenza = f'{zs} - (({nome_layer_dtm}@1)*1)'

                    # Definisci le voci del raster da calcolare
                    raster_dtm_entry = QgsRasterCalculatorEntry()
                    raster_dtm_entry.ref = f'{nome_layer_dtm}@1'
                    raster_dtm_entry.raster = raster_dtm
                    raster_dtm_entry.bandNumber = 1

                    # Esegui il calcolo della differenza
                    output_raster_differenza = f'{Path_cartella_lavoro}/Diff/Differenza_{counter}.tif'
                    calc = QgsRasterCalculator(espressione_differenza, output_raster_differenza,
                            'GTiff', raster_dtm.extent(), raster_dtm.width(),
                            raster_dtm.height(), [raster_dtm_entry])
                    calc.processCalculation()
                
                    # Carica il raster della differenza calcolato
                    ds = gdal.Open(output_raster_differenza)
                    band = ds.GetRasterBand(1)
                    data = band.ReadAsArray()

                    # Calcola la cella del pixel in metri quadrati (assumendo che il raster utilizzi coordinate metriche)
                    geotransform = ds.GetGeoTransform()
                    pixel_area = abs(geotransform[1] * geotransform[5])
                
                    # Calcola il volume sommando tutti i valori del raster differenza moltiplicati per l'area del pixel
                    volume = np.sum(data[(data > 0)& (data < zs)]) * pixel_area  # Considera solo i valori positivi per il volume
                    volume_formattato = float(format(volume, ".1f"))
                    # print(f"Volume = {volume_formattato} m^3")

                    # Chiudi il dataset gdal
                    ds = None
                    
                    volumi_invaso.append([zs, area_totale, volume_formattato])
                    layers_1.append(invaso_layer)
                    layer_da_rimuovere.append(invaso_layer)
                
                    
                    
                    #print(volume_formattato)
                    #print(volume_vasca)
                    errore = abs(volume_formattato-volume_vasca)
                    #print(errore)
                    #print(deltah)
                    #print(incremento)
                    if volume_formattato > (volume_vasca+10000) and counter2 == 1:
                        Vasca_remuve_name = f"Vasca_{counter}_{zs}"
                        # Ottieni il layer raster dal progetto
                        Vasca_layer_remuved = QgsProject.instance().mapLayersByName(Vasca_remuve_name)[0]
                        if Vasca_layer_remuved is not None:
                            # Ottieni il percorso del file raster
                            Vasca_remuved_path = Vasca_layer_remuved.dataProvider().dataSourceUri()
                            #Rimuovi layer dal progetto
                            QgsProject.instance().removeMapLayer(Vasca_layer_remuved.id())
                        break
                    elif volume_formattato < volume_vasca and counter3 < 2:
                        deltah += DHS
                        DHS = abs(DHS)
                        #print(f'V= {volume_formattato}')
                        #print(f'Vasca= {volume_vasca}')
                        #print(f'deltah< = {deltah}')
                        #print(DHS)
                    elif volume_formattato > volume_vasca:
                        #print(f'counter3 {counter3}')
                        deltah -= DHS/(counter3+1)
                        #print(f'deltah> = {deltah}')
                        #print(f'V= {volume_formattato}')
                        #print(f'Vasca= {volume_vasca}')
                        #DHS = -((deltah+Q)/zs)/2
                        #print(f'DHS{DHS}')
                        counter3 += 1
                    elif volume_formattato < volume_vasca and counter3 >=2:
                        #print(f'counter3 < {counter3}')
                        deltah += DHS/(counter3+1)
                        DHS = abs(DHS)
                        #print(f'V= {volume_formattato}')
                        #print(f'Vasca= {volume_vasca}')
                        #print(f'deltah< = {deltah}')
                        #print(DHS)
                        counter3 += 1
                    elif volume_formattato == volume_vasca:
                        break
                    print(zs)
                    #funziona bisogna farlo convergere al risultato il piu velocemente possibile 
                    counter2 += 1
                    Vasca_remuve_name = f"Vasca_{counter}_{zs}"
                    # Ottieni il layer raster dal progetto
                    Vasca_layer_remuved = QgsProject.instance().mapLayersByName(Vasca_remuve_name)[0]
                    if Vasca_layer_remuved is not None:
                        # Ottieni il percorso del file raster
                        Vasca_remuved_path = Vasca_layer_remuved.dataProvider().dataSourceUri()
                        #Rimuovi layer dal progetto
                        QgsProject.instance().removeMapLayer(Vasca_layer_remuved.id())
                
                if counter2 < ncicli and volume_min < volume_formattato < volume_max:
                    print(f'Ho trovato il mio sbrramento')
                elif counter2 == ncicli:
                    print(f'il ciclo non ha trovato una soluzione accettabile') #cambiare il numero di ciclci per permettere di andare oltre gli 8 metri
                print(f'volume calcolato = {volume_formattato}& zs = {zs}')
                volumi_filtrati = [riga for riga in volumi_invaso if volume_min < riga[2] < volume_max]
                
                #Salva le modifiche
                layer_sez_chiusura.commitChanges()
                
                #IL RIMUOVONE!! 
                remuve_name = f"Area_{counter}"
                # Ottieni il layer raster dal progetto
                layer_remuved = QgsProject.instance().mapLayersByName(remuve_name)[0]
                if layer_remuved is not None:
                    # Ottieni il percorso del file raster
                    remuved_path = layer_remuved.dataProvider().dataSourceUri()
                    #Rimuovi layer dal progetto
                    QgsProject.instance().removeMapLayer(layer_remuved.id())
        
                B_remuve_name = f"Bacino_{counter}"
                # Ottieni il layer raster dal progetto
                B_layer_remuved = QgsProject.instance().mapLayersByName(B_remuve_name)[0]
                if B_layer_remuved is not None:
                    # Ottieni il percorso del file raster
                    B_remuved_path = B_layer_remuved.dataProvider().dataSourceUri()
                    #Rimuovi layer dal progetto
                    QgsProject.instance().removeMapLayer(B_layer_remuved.id())
        
                V_remuve_name = f"Vect_{counter}"
                # Ottieni il layer raster dal progetto
                V_layer_remuved = QgsProject.instance().mapLayersByName(V_remuve_name)[0]
                if V_layer_remuved is not None:
                    # Ottieni il percorso del file raster
                    V_remuved_path = V_layer_remuved.dataProvider().dataSourceUri()
                    #Rimuovi layer dal progetto
                    QgsProject.instance().removeMapLayer(V_layer_remuved.id())
        
                    
                #RIMUOVO I LAYER CHE NON RIENTRANO CON LE CARATTERISTICHE DELL'INVASO
                
                project = QgsProject.instance()
                layer_prefisso = f"Invaso_{counter}_"
                valore_lista = volumi_filtrati
                if not valore_lista:
                    #se la lista è vuota cancella tutti i layer con il prefisso invaso_x_
                    #layer_da_rimuovere = [layer for layer in layers.values() if layer.name().startswith(layer_prefisso)]
                    for layer in layer_da_rimuovere:
                        project.removeMapLayer(layer)
                    print(f"Rimosso {len(layer_da_rimuovere)} layer con prefisso '{layer_prefisso}'.")
                else:
                    #se la lista è riempita, seleziona i layer da tenere e rimuovi gli altri
                    valori_da_selezionare = [str(value[0]) for value in valore_lista]
                    layer_da_rimuovere = []
                    
                    for layer in layers_1:#layers.values():
                        if layer.name().startswith(layer_prefisso):
                            suffisso = layer.name().replace(layer_prefisso, '')
                        
                            if suffisso not in valori_da_selezionare:
                                layer_da_rimuovere.append(layer)
                    for layer in layer_da_rimuovere:
                        project.removeMapLayer(layer)
                    print(f"Ho rimosso {len(layer_da_rimuovere)}. Layer non nella lista Invaso_... {valori_da_selezionare}")
                
                #salvataggio zs/area/volume
                if volume_formattato < volume_max:
                    layer_change0 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change0.dataProvider()
                    zs_bac_idx = layer_change0.fields().lookupField('zs')
                    zs_formatted = f"{zs:.4f}"
                    att = {zs_bac_idx:zs_formatted} 
                    feat = layer_change0.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                    layer_change1 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change1.dataProvider()
                    area_bac_idx = layer_change1.fields().lookupField('Area Bac')
                    att = {area_bac_idx: f'{area_totale}'}
                    feat = layer_change1.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                    layer_change2 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change2.dataProvider()
                    volume_bac_idx = layer_change2.fields().lookupField('Volume Bac')
                    att = {volume_bac_idx: f'{volume_formattato}'}
                    feat = layer_change2.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                elif volume_formattato > volume_max:
                    layer_change0 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change0.dataProvider()
                    zs_bac_idx = layer_change0.fields().lookupField('zs')
                    att = {zs_bac_idx: f"{Q}"}
                    feat = layer_change0.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                    layer_change1 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change1.dataProvider()
                    area_bac_idx = layer_change1.fields().lookupField('Area Bac')
                    att = {area_bac_idx: 0}
                    feat = layer_change1.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                    layer_change2 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change2.dataProvider()
                    volume_bac_idx = layer_change2.fields().lookupField('Volume Bac')
                    att = {volume_bac_idx: 0}
                    feat = layer_change2.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                elif volume_formattato < volume_min:
                    layer_change0 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change0.dataProvider()
                    zs_bac_idx = layer_change0.fields().lookupField('zs')
                    att = {zs_bac_idx: f"{Q}"}
                    feat = layer_change0.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                
                    layer_change1 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change1.dataProvider()
                    area_bac_idx = layer_change1.fields().lookupField('Area Bac')
                    att = {area_bac_idx: 0}
                    feat = layer_change1.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                    
                    layer_change2 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                    prov = layer_change2.dataProvider()
                    volume_bac_idx = layer_change2.fields().lookupField('Volume Bac')
                    att = {volume_bac_idx: 0}
                    feat = layer_change2.getFeature(counter-1)
                    prov.changeAttributeValues({feat.id():att})
                
                
                
                layer_change5 = QgsProject().instance().mapLayersByName('Sez_Chiusura')[0]
                prov = layer_change5.dataProvider()
                HS_bac_idx = layer_change5.fields().lookupField('HS')
                if volume_min < volume_formattato < volume_max:
                    hs = (zs - Q)
                    #print(f'HS{hs}')
                    att = {HS_bac_idx: f'{hs:.2f}'}
                else:
                    hs = 0
                    #print(f'HS{hs}')
                    att = {HS_bac_idx: f'{hs}'}
                feat = layer_change5.getFeature(counter-1)
                prov.changeAttributeValues({feat.id():att})
                
                DHS = DHSu #mi serve per resettare il ciclo while
                counter +=1
            
        else:
            print("Il layer non è di tipo punto.")
                
    
#Applica la funzione al layer
apply_upslope_area_to_points(layer_sez_chiusura, DHSu)

# valutazione sbarramento
layer_name = "Sez_Chiusura"
# Ottieni il layer punti esistente
layers = QgsProject.instance().mapLayersByName(layer_name)
if not layers:
    raise ValueError(f"Layer '{layer_name}' non trovato nel progetto corrente.")
layer = layers[0]
# Definisci la color ramp di rossi
color_ramp = QgsGradientColorRamp(QColor(10, 210, 10 ), QColor(255, 0, 0)) #dal rosso al verde
# Imposta il campo per la classificazione
field_name = "HS"
# Verifica se il campo esiste nel layer
if field_name not in [field.name() for field in layer.fields()]:
    raise ValueError(f"Campo '{field_name}' non trovato nel layer '{layer_name}'.")
# Calcola i valori min e max per la classificazione
layer_min = layer.minimumValue(layer.fields().lookupField(field_name))
layer_max = layer.maximumValue(layer.fields().lookupField(field_name))
# Numero di classi
num_classes = count #numero di punti
# Calcola gli intervalli per le classi
step = ((layer_max - layer_min) / num_classes)
ranges = []
for i in range(num_classes):
    lower_value = layer_min + (i * step)
    upper_value = layer_min + ((i + 1) * step)
    symbol = QgsSymbol.defaultSymbol(layer.geometryType())
    symbol.setColor(color_ramp.color(float(i) / num_classes))
    label = f"{lower_value:.2f} - {upper_value:.2f}"
    ranges.append(QgsRendererRange(lower_value, upper_value, symbol, label))

# Crea il renderer graduato
renderer = QgsGraduatedSymbolRenderer(field_name, ranges)
renderer.setMode(QgsGraduatedSymbolRenderer.EqualInterval)
# Applica il renderer al layer
layer.setRenderer(renderer)

# Assumi che il layer vettoriale sia il primo layer nella lista
layer_name = 'Sez_Chiusura'
layer = QgsProject.instance().mapLayersByName(layer_name)[0]

# Configura le impostazioni per l'etichetta
label_settings = QgsPalLayerSettings()
label_settings.isExpression = True
label_settings.fieldName = "HS"  # Sostituisci con il nome del campo che vuoi etichettare
label_settings.placement = QgsPalLayerSettings.OverPoint

# Configura il formato del testo
text_format = QgsTextFormat()
text_format.setFont(QFont("Arial", 12))
text_format.setSize(12)
text_format.setColor(QColor(0, 0, 0))  # Colore rosso per il testo

# Configura il buffer per le etichette
buffer_settings = QgsTextBufferSettings()
buffer_settings.setEnabled(True)
buffer_settings.setSize(1)  # Dimensione del buffer
buffer_settings.setColor(QColor(255, 255, 255))  # Colore del buffer (bianco)

# Applica il buffer al formato del testo
text_format.setBuffer(buffer_settings)

label_settings.setFormat(text_format)

# Imposta l'offset delle etichette
label_settings.xOffset = 5  # Spostamento a destra
label_settings.yOffset = -5  # Spostamento verso l'alto

# Applica le impostazioni dell'etichetta al layer
layer.setLabelsEnabled(True)
layer.setLabeling(QgsVectorLayerSimpleLabeling(label_settings))
layer.triggerRepaint()


print("############################")
print("########FINE CALCOLO########")
print("############################")