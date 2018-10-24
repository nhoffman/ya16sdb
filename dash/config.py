import multiprocessing

bind = "127.0.0.1:8005"  # 0.0.0.0:443  # Docker
certfile = 'test.crt'
# chdir = '/dash'  # docker
keyfile = 'test.key'
threads = 4
workers = 4  # multiprocessing.cpu_count() * 2 + 1
