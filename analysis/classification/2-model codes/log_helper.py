#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Simple logging helper: redirect stdout/stderr to a timestamped log file under results/logs

Usage: at the top of your script do:
    from log_helper import setup_logging
    setup_logging(__file__)

This will duplicate all print() output to both console and the log file.
"""
import os
import sys
from datetime import datetime

def setup_logging(script_path=None, prefix=None):
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_root = os.path.dirname(script_dir)
        logs_dir = os.path.join(project_root, 'results', 'logs')
        os.makedirs(logs_dir, exist_ok=True)

        if script_path:
            name = os.path.splitext(os.path.basename(script_path))[0]
        else:
            name = prefix or 'script'

        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        logfile = os.path.join(logs_dir, f"{name}_{timestamp}.log")

        class Tee:
            def __init__(self, *streams):
                self.streams = streams

            def write(self, data):
                for s in self.streams:
                    s.write(data)

            def flush(self):
                for s in self.streams:
                    try:
                        s.flush()
                    except Exception:
                        pass

        log_f = open(logfile, 'a', encoding='utf-8')
        sys.stdout = Tee(sys.__stdout__, log_f)
        sys.stderr = Tee(sys.__stderr__, log_f)

        print(f"[LOG] 输出同时记录到: {logfile}")
        return logfile
    except Exception as e:
        print(f"无法创建日志文件: {e}")
        return None
