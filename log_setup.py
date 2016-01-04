import logging
import logging.config

log_conf = {
    'version': 1,
    'filters': [],
    'formatters': {
        'simple': {
            'format': '%(levelname)s  %(message)s'
        },
        'detailed': {
            'class': 'logging.Formatter',
            'format': '%(levelname)-8s %(name)s:%(filename)s:%(lineno)d  %(message)s',
        },
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
        },
    },
    'handlers': {
        'console': {
            'class': 'logging.StreamHandler',
            'formatter': 'simple',
            'level': 'INFO',
        },
        'file': {
            'class': 'logging.FileHandler',
            'filename': 'processing.log',
            'mode': 'w',
            'level': 'DEBUG',
            'formatter': 'detailed',
        },
        'errors': {
            'class': 'logging.FileHandler',
            'filename': 'errors.log',
            'mode': 'w',
            'level': 'ERROR',
            'formatter': 'detailed',
        },
    },
    'loggers': {
        'foo': {
            'handlers': ['file']
        },

        'root': {
            'level': 'DEBUG',
            'handlers': ['console', 'file', 'errors'],

        },
    },
}


def getLogger(name='root'):
    logging.config.dictConfig(log_conf)
    logger = logging.getLogger(name)
    return logger


