#!/usr/bin/env python3

import argparse
import re
from collections import namedtuple
from itertools import groupby
import os
import html
from fnmatch import fnmatch

from lxml import etree

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

def get_clang_tidy_warnings(input):
    files, items = parse_clang_tidy_output(input)

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
    for f in ref_files:
        items_by_file[f] = []

    for item in items:
        items_by_file[os.path.normpath(item.file)].append(item)

    def get_comment(file, lineno, cr):

        items = [i for i in items_by_file[os.path.normpath(file)] if i.lineno == lineno]
        if len(items) == 0:
            return None

        msgs = []

        for item in items:
            msg = html.escape(item.itemtype + ": "+item.msg+"\n"+item.code)
            msg = "\n".join(['<pre style="display:block;white-space:pre-wrap;">{}</pre>'.format(l) for l in msg.split("\n")])
            msg = add_filelinks(msg, cr)
            msgs.append(msg)

        return msgs

    # ath_prefix = os.path.commonprefix(files)
    # print(path_prefix)

    return files, get_comment

class CppcheckItem:
    def __init__(self, file, line, verbose, severity, id):
        self.file = file
        self.line = line
        self.verbose = verbose
        self.severity = severity
        self.id = id

    def __str__(self):
        return "{severity}: {msg} [{id}]".format(
            severity= self.severity,
            msg= self.verbose,
            id= self.id
        )

def get_cppcheck_warnings(input):
    root = etree.fromstring(input)

    errors = root.xpath("/results/errors")[0]

    items_by_file = {}


    for error in errors:
        location = error.xpath("./location")[0]
        item = CppcheckItem(
            file = location.attrib["file"],
            line = int(location.attrib["line"]),
            verbose = error.attrib["verbose"],
            severity = error.attrib["severity"],
            id = error.attrib["id"]
        )

        if not item.file in items_by_file:
            items_by_file[item.file] = {}

        if not item.line in items_by_file[item.file]:
            items_by_file[item.file][item.line] = []
        
        items_by_file[item.file][item.line].append(item)

    def get_comment(file, lineno, cr):
        if file in items_by_file and lineno in items_by_file[file]:
            items = items_by_file[file][lineno]
        else:
            return None

        fmt = '<pre style="white-space:pre-wrap;display:block;">{}</pre>'
        return list(map(lambda item: fmt.format(str(item)), items))

    return list(items_by_file.keys()), get_comment



def main():
    p = argparse.ArgumentParser()
    p.add_argument("mode", choices=("clang-tidy", "cppcheck"))
    p.add_argument("inputfile")
    p.add_argument("reportdir", default="report")
    p.add_argument("--exclude", "-e", action="append", default=[])

    args = p.parse_args()


    if args.mode == "clang-tidy":
        with open(args.inputfile, "r", encoding="utf-8") as f:
            inputstr = f.read()
        files, get_comment = get_clang_tidy_warnings(inputstr)
        title = "ACTS clang-tidy report"
    elif args.mode == "cppcheck":
        with open(args.inputfile, "rb") as f:
            inputstr = f.read()
        files, get_comment = get_cppcheck_warnings(inputstr)
        title = "ACTS cppcheck report"

    files = filter(lambda f: not any(fnmatch(f, e) for e in args.exclude), files)

    cr = CodeReport(files, 
                    title=title,
                    get_comment=get_comment)


    for file, content in cr:
        with open(os.path.join(args.reportdir, file), "w", encoding='utf-8') as f:
            print(os.path.join(args.reportdir, file))
            f.write(content)




if "__main__" == __name__:
    main()
