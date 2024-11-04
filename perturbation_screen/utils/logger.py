import logging


def setup_logger(name: str, log_file: str, level=logging.INFO):
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )

    handler = logging.FileHandler(log_file)
    handler.setFormatter(formatter)

    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger


def log(msg: str, logger: logging.Logger | None = None, level: int = logging.INFO):
    if logger:
        logger.log(level, msg)
