import hdf5_to_dict as io
import sys

preferred_keys = ['N1','N2','N3',
                  'a',
                  'M_unit','B_unit','L_unit','Ne_unit']

if len(sys.argv) != 2:
  print('ERROR format is')
  print('  python hdr.py [filename]')
  sys.exit()

fnam = sys.argv[1]
hdr = io.load_hdr(fnam)

used_keys = []
def print_val(vnam):
  if vnam in hdr.keys() and vnam not in used_keys:
    format_string = '%e' if type(hdr[vnam]) is float else '%s'
    print(vnam,'=',format_string % hdr[vnam])
    used_keys.append(vnam)

for k in preferred_keys:
  print_val(k)

for k in hdr.keys():
  print_val(k)
