from time import localtime, strftime

def log_message(*args, **kwargs):
    print(
        '[' + strftime('%H:%M:%S', localtime()) + ']',
        *args, **kwargs)
