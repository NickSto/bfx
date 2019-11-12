__version__ = '0.2'
from collections import OrderedDict

FLAGS = OrderedDict((
  (1, 'read paired'),
  (2, 'read mapped in proper pair'),
  (4, 'read unmapped'),
  (8, 'mate unmapped'),
  (16, 'read reverse strand'),
  (32, 'mate reverse strand'),
  (64, 'first in pair'),
  (128, 'second in pair'),
  (256, 'not primary alignment'),
  (512, 'read fails platform/vendor quality checks'),
  (1024, 'read is PCR or optical duplicate'),
  (2048, 'supplementary alignment'),
))

def decompose(flags_int):
  flag_statuses = {}
  for flag in FLAGS:
    if flags_int & flag:
      flag_statuses[flag] = True
    else:
      flag_statuses[flag] = False
  return flag_statuses