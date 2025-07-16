#!/usr/bin/env python

"""Script to check for missing function, class and class method documentation in the
source/api documentation files.

Author: Samuel Rowlinson
"""

import os
from collections import OrderedDict
import finesse


class _parsed_docs:
    def __init__(self):
        self.classes = []
        self.class_methods = {}
        self.functions = []


def iterate_class_methods(file, classes):
    classes_w_methods = {}
    class_docs = []
    for i, (spos, name, is_abstract, is_private) in enumerate(classes):
        if is_private:
            continue
        file.seek(spos)
        classes_w_methods[name] = set()
        if spos == classes[-1][0]:
            allclass = file.read()
        else:
            allclass = file.read(classes[i + 1][0] - spos)
        classes_splitlines = list(iter(allclass.splitlines()))
        for j, cline in enumerate(classes_splitlines):
            stripped = cline.lstrip()
            if not j:
                class_docs.append(
                    (name, (stripped.startswith('"""') or stripped.startswith('r"""')))
                )
            # only grab methods within one level of indentation (class methods)
            if stripped.startswith("def") and len(cline) - len(stripped) == 4:
                meth_name = (stripped.split(" ")[1]).split("(")[0]
                is_property = classes_splitlines[j - 1].lstrip().startswith("@property")
                # grab only public methods - don't add __init__ if class is abstract
                if not meth_name.startswith("_") or (
                    meth_name == "__init__" and not is_abstract
                ):
                    # does the method have docstrings
                    if classes_splitlines[j + 1].lstrip().startswith(
                        '"""'
                    ) or classes_splitlines[j + 1].lstrip().startswith('r"""'):
                        classes_w_methods[name].add((meth_name, True, is_property))
                    else:
                        # protect against class properties being added twice
                        if not (meth_name, True, True) in classes_w_methods[name]:
                            classes_w_methods[name].add((meth_name, False, is_property))
    return class_docs, classes_w_methods


def iterate_functions(file, functions):
    funcs = []
    for spos, name in functions:
        file.seek(spos)
        docline = file.readline()
        if not name.startswith("_"):
            funcs.append(
                (
                    name,
                    (
                        docline.lstrip().startswith('"""')
                        or docline.lstrip().startswith('r"""')
                    ),
                )
            )
    return funcs


def classes_methods_in(file):
    classes = []  # holds tuples of (streampos, class_name, is_abstract, is_private)
    functions = []  # holds tuples of (streampos, function_name)
    with open(file) as f:
        line = f.readline()
        in_mod_docstrings = False
        while line:
            if line.startswith('"""') and not in_mod_docstrings:
                in_mod_docstrings = True
            elif line.startswith('"""') and in_mod_docstrings:
                in_mod_docstrings = False
            # don't consider any "class" string found in docstrings
            if line.startswith("class") and not in_mod_docstrings:
                name_base = (line.split(" ")[1]).split("(")
                classes.append(
                    (
                        f.tell(),
                        name_base[0].split(":")[0],
                        "ABC" in name_base[1] if len(name_base) > 1 else False,
                        name_base[0].startswith("_"),
                    )
                )
            # check for free function docstrings
            elif line.startswith("def") and not in_mod_docstrings:
                name = (line.split(" ")[1]).split("(")[0]
                functions.append((f.tell(), name))
            line = f.readline()
        pd = _parsed_docs()
        pd.classes, pd.class_methods = iterate_class_methods(f, classes)
        pd.functions = iterate_functions(f, functions)
    return pd


def undocumented_functions(functions, docfile_contents, undoc_items):
    for name, has_docstrings in functions:
        # (type, is function included in docfile, does function have docstrings)
        if name not in docfile_contents:
            undoc_items[name] = ("Function", False, has_docstrings)
        else:
            if not has_docstrings:
                undoc_items[name] = ("Function", True, False)


def undocumented_classes(classes, docfile_contents, undoc_items):
    for name, has_docstrings in classes:
        if name not in docfile_contents:
            undoc_items[name] = ("Class", False, has_docstrings)
        else:
            if not has_docstrings:
                undoc_items[name] = ("Class", True, False)


def undocumented_class_methods(class_methods, docfile_contents, undoc_items):
    for class_name, methods in class_methods.items():
        for name, has_docstrings, is_property in methods:
            pre = "Property" if is_property else "Method"
            if not class_name + "." + name in docfile_contents:
                undoc_items[class_name + "." + name] = (pre, False, has_docstrings)
            else:
                if not has_docstrings:
                    undoc_items[class_name + "." + name] = (pre, True, False)


def check_doc_file(file, parsed_docs):
    if not os.path.isfile(file):
        print("\x1b[0;30;41mNo documentation file: {} present\x1b[0m".format(file))
        return False
    with open(file) as f:
        docfile_contents = f.read()
        # stores [name, (type_str, name in docfile, name has docstrings)]
        undoc_items = OrderedDict()
        undocumented_functions(parsed_docs.functions, docfile_contents, undoc_items)
        undocumented_classes(parsed_docs.classes, docfile_contents, undoc_items)
        undocumented_class_methods(
            parsed_docs.class_methods, docfile_contents, undoc_items
        )
        for item, (pre, is_in_docfile, has_docstrings) in undoc_items.items():
            if not is_in_docfile:
                if not has_docstrings:
                    print(
                        f"\x1b[0;31;40m{pre}: {item} not documented in file: {f.name}\x1b[0m"
                    )
                else:
                    print(f"{pre}: {item} not documented in file: {f.name}")
            else:
                print(
                    f"\x1b[0;33;40m{pre}: {item} included in documentation but "
                    "has no docstrings\x1b[0m"
                )
    return len(undoc_items)


def check_doc_files(files, parsed_docs, src_to_doc_path):
    all_undoc_items = []
    for file in files:
        if not os.path.isfile(file):
            print("\x1b[0;30;41mNo documentation file: {} present\x1b[0m".format(file))
            continue
        with open(file) as f:
            docfile_contents = f.read()
            # stores [name, (type_str, name in docfile, name has docstrings)]
            undoc_items = OrderedDict()
            undocumented_functions(parsed_docs.functions, docfile_contents, undoc_items)
            undocumented_classes(parsed_docs.classes, docfile_contents, undoc_items)
            undocumented_class_methods(
                parsed_docs.class_methods, docfile_contents, undoc_items
            )
            all_undoc_items.append(undoc_items)
    if not all_undoc_items:
        return 0
    undoc_items_final = OrderedDict()
    for item in all_undoc_items[0].keys():
        item_in_all = True
        for j in range(1, len(all_undoc_items)):
            if item not in all_undoc_items[j].keys() or all_undoc_items[j][item][1]:
                item_in_all = False
                break
        if item_in_all:
            undoc_items_final[item] = all_undoc_items[0][item]
    for item, (pre, is_in_docfile, has_docstrings) in undoc_items_final.items():
        if not is_in_docfile:
            if not has_docstrings:
                print(
                    f"\x1b[0;31;40m{pre}: {item} not documented in any files of "
                    f"directory: {src_to_doc_path}\x1b[0m"
                )
            else:
                print(
                    f"{pre}: {item} not documented in any files of directory: {src_to_doc_path}"
                )
        else:
            print(
                f"\x1b[0;33;40m{pre}: {item} included in documentation but "
                "has no docstrings\x1b[0m"
            )
    return len(undoc_items_final)


def check_single(src_file, doc_file):
    pd = classes_methods_in(src_file)
    if check_doc_file(doc_file, pd):
        print()


def check_multiple(src_file, doc_files, src_to_doc_path):
    pd = classes_methods_in(src_file)
    if check_doc_files(doc_files, pd, src_to_doc_path):
        print()


def recurse_check(file, pre_path=""):
    if file.endswith(".py") and file != "__init__.py":
        src_file = os.path.abspath(
            finesse.__file__.replace("__init__.py", pre_path + file)
        )
        doc_real_path = os.path.dirname(os.path.realpath(__file__))
        src_to_doc_path = (
            doc_real_path
            + f"/source/api/{pre_path + '/' if pre_path else ''}{file.split('.')[0]}"
        )
        this_doc_file = doc_real_path + "/source/api{0}finesse.{1}{2}.rst".format(
            "/" + pre_path if pre_path else "/",
            pre_path.split("/")[0] + "." if pre_path else "",
            file.split(".")[0],
        )
        # check whether the documentation for the corresponding
        # source file is split into its own sub-directory
        if os.path.isdir(src_to_doc_path):
            sub_doc_files = [this_doc_file]
            for sub_file in os.listdir(src_to_doc_path):
                if sub_file.endswith(".rst"):
                    sub_doc_files.append(src_to_doc_path + "/" + sub_file)
            check_multiple(src_file, sub_doc_files, src_to_doc_path)
        else:
            check_single(src_file, this_doc_file)
    elif "." not in file and not file.startswith("_"):
        for sub_file in os.listdir(
            finesse.__file__.replace("__init__.py", pre_path + file + "/")
        ):
            recurse_check(sub_file, pre_path + file + "/")


def check_all():
    for file in os.listdir(finesse.__file__.replace("__init__.py", "")):
        recurse_check(file)


print("Key")
print("---")
print(
    "\tLines highlighted \x1b[0;30;41min this style\x1b[0m indicate a non-existent documentation file."
)
print(
    "\tClass-methods highlighted in \x1b[0;31;40mred\x1b[0m do not have docstrings and are "
    "not included in their associated documentation file."
)
print(
    "\tClass-methods highlighted in \x1b[0;33;40myellow\x1b[0m do not have docstrings but are "
    "included in the documentation."
)
print(
    "\tClass-methods highlighted in white have docstrings but are not "
    "included in their associated documentation file.\n"
)
print("-------------------------------------------------\n")
check_all()
