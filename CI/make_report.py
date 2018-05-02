#!/usr/bin/env python3

import argparse
import re
from collections import namedtuple
from itertools import groupby
import os
import html

from codereport import CodeReport


class Item:
    def __init__(self, file, lineno, col, itemtype, msg, code):
        self.file = file
        self.lineno = lineno
        self.col = col
        self.itemtype = itemtype
        self.msg = msg
        self.code = code

    def __eq__(self, other):
        return (self.file == other.file 
               and self.lineno == other.lineno 
               and self.col == other.col)

    def __hash__(self):
        return hash((self.file, self.lineno, self.col))


def parse_clang_tidy_item(itemstr):

    m = re.match(r"([/.\-\w]+):(\d+):(\d+): (.*?):(.*?)\n((?:.|\n)*)", itemstr)

    if m is None:
        print(itemstr)

    item = Item(
        file = m.group(1),
        lineno = int(m.group(2)),
        col = int(m.group(3)),
        itemtype = m.group(4),
        msg = m.group(5),
        code = m.group(6),
    )

    return item



def parse_clang_tidy_output(output):

    # cleanup
    itemstr = re.sub(r"Enabled checks:\n(?:.|\n)+\n\n", "", output)
    itemstr = re.sub(r"clang-tidy-\d\.\d.*\n?", "", itemstr)
    itemstr = re.sub(r".*-header-filter.*", "", itemstr)

    items = []
    prevstart = 0

    matches = list(re.finditer(r"([\w/.\-]+):(\d+):(\d+): warning:", itemstr))
    for idx, m in enumerate(matches):
        # print(m)
        start, end = m.span()
        if idx > 0:
            item = itemstr[prevstart:start]
            items.append(item)
        if idx+1 == len(matches):
            item = itemstr[start:]
            items.append(item)
        prevstart = start

    items = set(map(parse_clang_tidy_item, items))

    files = set(map(lambda item: item.file, items))

    return files, items


def main():
    p = argparse.ArgumentParser()
    p.add_argument("inputfile", type=argparse.FileType("r"))
    p.add_argument("reportdir", default="report")

    args = p.parse_args()


    files, items = parse_clang_tidy_output(args.inputfile.read())

    files = filter(os.path.exists, files)
    files = list(files)
    items = list(filter(lambda i: os.path.exists(i.file), items))

    ref_files = list(map(os.path.normpath, files))

    def add_filelinks(msg, cr):
        def repl(m):
            if os.path.normpath(m.group(1)) in ref_files:
                href = cr.get_file_link(m.group(1), m.group(2), m.group(3))
                return '<a href="{}">{}</a>'.format(href, m.group(0))
            else:
                return m.group(0)

        return re.sub(r"([.\/\w]+):(\d+):(\d+):", repl, msg)

    items_by_file = {}
    for f in files:
        items_by_file[f] = []

    for item in items:
        items_by_file[item.file].append(item)

    def get_comment(file, lineno, cr):
        items = [i for i in items_by_file[file] if i.lineno == lineno]
        if len(items) == 0:
            return None

        msgs = []

        for item in items:
            msg = html.escape(item.itemtype + ": "+item.msg+"\n"+item.code)
            msg = "\n".join(['<pre style="display:block;white-space:pre-wrap;">{}</pre>'.format(l) for l in msg.split("\n")])
            msg = add_filelinks(msg, cr)
            msgs.append(msg)

        return msgs

    path_prefix = os.path.commonprefix(files)
    print(path_prefix)

    cr = CodeReport(files, 
                    title="ACTS clang-tidy report",
                    get_comment=get_comment)


    for file, content in cr:
        with open(os.path.join(args.reportdir, file), "w", encoding='utf-8') as f:
            print(os.path.join(args.reportdir, file))
            f.write(content)




if "__main__" == __name__:
    main()
