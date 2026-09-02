
import ijson
path = r'data\clean\species\chenopodium_quinoa_all_sources.json'
with open(path, 'rb') as f:
    found = False
    for prefix, event, value in ijson.parse(f):
        if event == 'map_key' and value == 'genes' and prefix == '':
            found = True
            continue
        if found:
            print('Premier événement après genes:', event)
            break