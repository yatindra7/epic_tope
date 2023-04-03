from aaindex import aaindex1
import json

with open('data/aaindex1.json', 'w') as aa_json:
    json.dump(aaindex1.parse_aaindex(), aa_json)