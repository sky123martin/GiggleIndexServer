import subprocess
from config import config
import sqlite3
from sys import argv

config = config
conn = sqlite3.connect(config.DB)

# https://github.com/ryanlayer/cmd2web/blob/master/src/server.py line 173
#    out_file_name = '/tmp/' + str(random.randint(0,sys.maxsize)) + '.out'
#     f = open(out_file_name, 'w')
#     try:
#         proc = subprocess.check_call(cmd,
#                                      stderr=None,
#                                      stdout=f,
#                                      timeout=timeout)
#         sys.stderr.write("\n\n\Command: {0}\n\n\n".format(cmd))
#         # res = subprocess.check_output(cmd,stderr=sys.stderr)
#         # sys.stderr.write("\n\n\nResult: {0}\n\n\n".format(res))
#     except subprocess.TimeoutExpired as e:
#         print('Time Out')
#         logger.error("Time limit for current request exceed.")
#         return cmd2web.Server.error('Time limit for current request exceed.')
#     except Exception as e:
#         return cmd2web.Server.error(str(e))

#     f.close()
