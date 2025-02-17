# Cómo usar `blastoise.py`

## Introducción

Al momento de escribir este documento (04/12/2024) el software, aunque funcional, está en la etapa final, en donde tiene que ser preparado a nivel de documentación como de disponibilidad. Con esto me refiero a que, aunque pueda ser utilizado correctamente, no tiene preparada la correcta documentación ni en GitHub, ni en Docker Hub, etc.

Por ello espero poder explicar su funcionamiento y utilización de forma correcta.

## Requerimientos

Para poder utilizar el programa necesitarás:

- Linux OS Ubuntu (lo he probado tanto en la versión 22.04 como la 24.04).
- Conda (herramienta de organización de environments). Yo utilizo la versión 24.5.0, pero supongo que otra versión valdrá así como Miniconda o Micromamba.

Hay otro método de uso, que es utilizando construcción de imágenes mediante **Docker Engine**, pero quizás pueda resultar más complicado y todavía no tengo tan trabajada su preparación para facilidad de uso.

## 1º paso: preparación del environment.

Dentro del repositorio, toda la información del programa está dentro de la carpeta [blastoise](../blastoise), la otra carpeta [docs](../docs) está dedicada solo a documentación.

En esa carpeta nos interesa 2 ficheros:
- [environment.yml](../blastoise/environment.yml) --> lo utilizaremos para crear el _environment_ con conda.
- [blastoise.py](../blastoise/blastoise.py) --> El programa principal. Todavía está sujeto a ser cambiado de nombre.

El primer paso será tener todos los paquetes y requerimientos del programa, para lo cual usaremos "conda" en la terminal de _bash_:

1. En la terminal de *bash* nos movemos dentro de "path" en donde se encuentra el archivo [environment.yml](../blastoise/environment.yml), utilizando bash realizamos: `conda env create -f environment.yml`. Con esto crearemos un _environment_ con todas las dependencias necesarias, que se llamará blastoise (el nombre lo tiene, porque así está descrito dentro del archivo [environment.yml](../blastoise/environment.yml)). La instalación de las dependencias tardará un rato en instalarse.


2. Una vez terminado, activamos nuestro nuevo _environment_ con `conda activate blastoise`. A partir de aquí ya tenemos las herramientas necesarias para ejecutar el programa.

## 2º paso: ejecución de blastoise.py

A continuación ejecutaremos el programa utilizando python, el cual se habrá instalado dentro de nuestro _environment_. Para ejecutarlo necesitamos dos archivos:

- El genoma a buscar en formato FASTA. Por ejemplo podría ser _Leishmania_infantim_whole_genome.fasta_.
- La secuencia de búsqueda en formato FASTA. Podría ser por ejemplo _ingi_hallmark.fasta_

Y llamamos al programa con esos dos elementos en _bash_. Recordamos que estamos en el "path" del archivo [blastoise.py](../blastoise/blastoise.py):

```bash
python blastoise.py \
-d <archivo secuencia de búsqueda FASTA> \
-g <archivo genoma FASTA>
```
Los "\\" son para saltar a la siguiente línea en _bash_ y así tener el código más ordenado.

Ejecutamos el código y nos preguntará por el nombre de la carpeta donde guardar los datos y, por el "path" en donde queremos guardar esa carpeta. Luego preguntará los siguientes parámetros:

- _Enter the identity for the first BLASTn step_:
  - Utilizad un 60%. Es la identidad que utilizará las búsquedas de BLASTn repetidas, usaremos este número basándonos en bibliografía tanto del equipo de Bringaud como del de Requena.
- _Enter the `word_size` value_:
  - Utilizad un 15. BLASTn por defecto usa 28 y en caso de Trypanosoma se usa 11, pero no hemos encontrado una diferencia significativa en usar 15 en este _software_, aunque sí que aumenta la velocidad de búsqueda.
- _Enter the `min_length` value_: 
  - Utilizad 100. Que es el valor mínimo de longitud para dar por válida una sequencia.
- _Enter the `extend_number` value_:
  - Utilizad 1000. Que es el valor de longitud de una secuencia que no llegue a los 1000 nt de longitud. Es decir, si mide 600 nt, será extendido por ambos lados hasta alcanzar los 1000.
- _Enter the number of the first run_:
  - Escribir 1. El valor es para enumerar las iteraciones del programa, en este caso empezará por "1" porque hemos escrito "1", si hubiéramos escrito "10", la primera iteración sería "10". Este valor es una reliquia del proceso de DEBUG del programa.

Una vez terminado, el programa se ejecutará, aquí depende mucho del HARDWARE que se tenga tanto de los procesos del ordenador que ejecute de entre medias. En mi caso con un **i9 de 11th generación y 64 GB de RAM** en el genoma de _L. infantum_ con el Hallmark de ingi, tardó casi 3h recién encendido, y en otra ocasión con los mismos hardwares, pero con el ordenador encendido por semanas, y por tanto, con múltiples procesos desconocidos abiertos, tardó 11 horas.

Una vez terminado tendréis todos los elementos de búsqueda en la carpeta que hayáis elegido. Dentro de esa carpeta tendréis la carpeta "execution_data" y dentro de esa carpeta habrá una archivo que se llama "Last_Data.csv" en donde tendrá todas las secuencias encontradas.

## 3º paso: depuración.

(05/12/24 en construcción)

El programa encuentra de forma iterativa todas las secuencias repetidas que buscamos, pero como tiene una habilidad de extensión, puede acabar encontrando otras secuencias y patrones que den lugar a otros elementos. Por lo tanto, al final nos encontramos con una amalgama de elementos repetidos en el genoma que son elementos SIDER y elementos repetidos no SIDER. Para poder filtrar bien los elementos tenemos los siguientes Scripts:

### 3.1 Union de datos en ambas hebras

[join_strands.py](../blastoise/extra/joinStrands.py)

"Last_Data.csv" tiene todas las secuencias repartidas tanto en el sentido positivo como negativo de las hebras. Esto es una repercusión de la forma que tiene BLASTn de encontrar las secuencias. Sin embargo, no podemos darle ese sentido debido a la estructura que tiene el genoma de _Leishmania_ (i.e., Directional Gene Clusters). Por ello, el primer paso es unir todos los elementos en el sentido "positivo", para ello usaremos [join_strands.py](../blastoise/extra/joinStrands.py):
  
  ```bash
  python ./extra/compare.py \
  <archivo CSV> \
  <archivo genoma FASTA> \
  plus
  ```
De aquí obtendremos un archivo denominado **merged_sequences.csv** como resultado.


### 3.2 Aplicación del test de SIDER

[true_sider_filter.py](../blastoise/extra/true_sider_filter.py)


Este _script_ se encargará de encontrar dentro del amalgama de elementos repetidos aquellos considerados realmente SIDER. Para ser considerado como SIDER en una búsqueda de cada secuencia en BLASTn frente al genoma original, consideramos que deben cumplir las siguientes premisas:
  - Aparecer en al menos 5 cromosomas diferentes con un valor esperado < 1.0E-09. Es decir, si $x$ es todas las secuencias de elementos SIDER, $C(x)$ es existir en al menos 5 cromosomas y $E(x)$ tener un valor esperado < 1.0E-09, entonces:

$$
\forall x \ (C(x) \ \land \ E(x) )
$$

Ejecutamos el script:

```bash
python ./extra/true_sider_filter.py \
-f <merged_secuences.csv> \
-d <archivo genoma FASTA> \
--word-size 15\
--recaught_file <archivo secuencia de búsqueda FASTA>
```

Una vez terminado, obtendremos dos archivos CSV:
  - **final_yes_data.csv** --> con aquellos elementos que han pasado el filtro. Este es el que nos interesa.
  - **final_no_data.csv** --> con aquellos elementos que no han pasado el filtro.

### 3.3. Corrección coordenadas

#### 3.3.1 Delimitación y fragmentación

[correct_coor_to_json.py](../blastoise/extra/correct_coor_to_json.py)

Este _script_ se encargará de mirar secuencia a secuencia y limitar las coordenadas de los elementos. El archivo resultante será un archivo tipo JSON con las secuencias y coordenadas resultantes.

```bash
python ./extra/correct_coor_to_json.py \
-f <path final_yes_data.csv> \
-d <archivo genoma FASTA> \
-o <path a la creación de la carpeta de destino deseada>
```

Dentro de la carpeta elegida, tendremos un archivo denominado **main_dict.json** con la información que necesitamos.

#### 3.3.2 True SIDER Test

[correct_coor_json_true_sider_test.py](../blastoise/extra/correct_coor_json_true_sider_test.py)

En este apartado, el _script _ utilizará los datos del archivo JSON **main_dict.json** para efectuar un análisis de True SIDER Test. La diferencia es que en vez de ejecutarse en el CSV, se utilizará el archivo JSON directamente.

```bash
python ./extra/correct_coor_json_true_sider_test.py \
-f <path main_dict.json> \
-o <path a la creación de la carpeta de destino deseada> \
-db <archivo genoma FASTA> 
```

Al finalizar el proceso, dentro de nuestra carpeta deseada tendremos un archivo JSON denominado **filtered_data.json** con la información necesaria.


#### 3.3.3 Obtención de secuencias

[correct_coor_json_to_csv.py](../blastoise/extra/correct_coor_json_to_csv.py)

Este _script_ realizará la transformación de los archivos JSON a un archivo CSV, de un uso más común.

```bash
python ./extra/correct_coor_json_to_csv.py \
-f  <path filtered_data.json> \
-o <path a la creación de la carpeta de destino deseada> \
-db <archivo genoma FASTA>\
-ndb <path final_no_data.csv>
```

Al finalizar el resultado tendremos dos archivos:
- **positive_database.csv** --> con los elementos deseados.
- **negative_database.csv** --> con los elementos rechazados, por si se necesitan estudiar en el futuro.