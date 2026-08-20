import bmrb
import json

if __name__ == '__main__':
  entry = bmrb.entry.fromFile('NEF_draft_2014_v6.star')
  entry.printTree()
  print json.dumps(entry, cls=bmrb.JSON_Encoder, indent=4, sort_keys=True)
