Aplikace pro tvorbu testovacích dat pro program DMU

vstupem programu jsou dva body ve WGS84 (EPSG:4326) souřadnicích a výstupem je složka s vygenerovanámi
testovacími daty a souborem pro program DMU.

## Instalace

```
git clone https://github.com/SteveeH/DMU_vstup_generator.git

cd DMU_vstup_generator

virtualenv venv
\venv\Scripts\activate

pip install -r requirements.txt
```

## Nastavení

nastavení parametrů je v souboru main.py dictionary settings např.

parametry:

- spacing - příčný a podélný rozestup bodů
- lines_count - počet podélných linií na obě strany od osy

```python
    settings = {
        "out_folder": "colas_vysocina",
        "axis_file_path": "axis.txt",
        "diff_file_path": "diff.txt",
        "dmt_file_path": "dmt.txt",
        "design_file_path": "desing.txt",
        "plg_file_path": "plg.txt",
        "kml_file_path": "output.kml",
        "start": Point(49.39032343930921, 15.604040568013803, default_print="wgs"),
        "end": Point(49.39067309938579, 15.603539879798616, default_print="wgs"),
        "spacing": 0.1,
        "lines_count": 40,
    }
```

## Použití

```
python main.py
```
