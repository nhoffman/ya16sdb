import multiprocessing

bind = "127.0.0.1:8005"  # 0.0.0.0:443  # Docker
certfile = 'test.crt'
# chdir = '/dash'  # docker
feather_file = '/mnt/disk15/molmicro/working/crosenth/src/ya16sdb/output/20181015/dedup/1200bp/named/filtered/filter_details.feather'
# feather_file = 'filter_details.feather'  # docker
keyfile = 'test.key'
threads = 4
workers = 4  # multiprocessing.cpu_count() * 2 + 1
