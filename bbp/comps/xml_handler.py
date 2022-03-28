#!/usr/bin/env python
"""
BSD 3-Clause License

Copyright (c) 2021, University of Southern California
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

This module implements functions that allows the Broadband Platform to
read and write workflows from/to a file.
"""
from __future__ import division, print_function

# Import Python modules
import os
import ast
import xml.dom.minidom as xml

# Import Broadband modules
import install_cfg
import validation_cfg
from module import Module

class Workflow(object):

    def __init__(self, module_list, val_obj, station_file, version):
        self.workflow = module_list
        self.val_obj = val_obj
        self.station_file = station_file
        self.version = version

'''
<BBP_Run_Specification>
        <Validation_Run event="Northridge" input_station_file="user/station_file.stl"/>
        <BBP_Modules>
                <BBP_Module>
                        <name>jbsim</name>
                        <staged_files>
                                <file>/path/to/input/file</file>
                                <file>/path/to/another/input/file</file>
                        </staged_files>
                        <arguments>
                                <argument type="string">foo</argument>
                                <argument type="int">123</argument>
                        </arguments>
                        <keyword_arguments>
                                <keyword_argument keyword="foo" type="string">bar</keyword_argument>
                </BBP_Module>
                <BBP_Module>
                        ...
                </BBP_Module>
        </BBP_Modules>
</BBP_Run_Specification>
'''

def expand_variables(argument):
    """
    This function expands variables, like $BBP_INSTALL, to create full
    path names
    """
    install = install_cfg.InstallCfg.getInstance()

    if argument.find('$BBP_INSTALL/') == 0:
        argument = argument.replace('$BBP_INSTALL/',
                                    "%s/" % install.A_INSTALL_ROOT)
    if argument.find('$BBP_INSTALL_GF/') == 0:
        argument = argument.replace('$BBP_INSTALL_GF/',
                                    "%s/" % install.A_GF_DIR)
    if argument.find('$BBP_INSTALL_VAL/') == 0:
        argument = argument.replace('$BBP_INSTALL_VAL/',
                                    "%s/" % install.A_VAL_DIR)
    if argument.find('$BBP_DATA_DIR/') == 0:
        argument = argument.replace('$BBP_DATA_DIR/',
                                    "%s/" % install.A_DATA_ROOT)

    return argument

def insert_variables(argument):
    """
    This function substitutes parts of a path with the corresponding
    BBP variable, like $BBP_INSTALL, if it finds a match
    """
    install = install_cfg.InstallCfg.getInstance()

    if argument.find("%s/" % install.A_INSTALL_ROOT) == 0:
        argument = argument.replace("%s/" % install.A_INSTALL_ROOT,
                                    '$BBP_INSTALL/')
    if argument.find("%s/" % install.A_GF_DIR) == 0:
        argument = argument.replace("%s/" % install.A_GF_DIR,
                                    '$BBP_INSTALL_GF/')
    if argument.find("%s/" % install.A_VAL_DIR) == 0:
        argument = argument.replace("%s/" % install.A_VAL_DIR,
                                    '$BBP_INSTALL_VAL/')
    if argument.find("%s/" % install.A_DATA_ROOT) == 0:
        argument = argument.replace("%s/" % install.A_DATA_ROOT,
                                    '$BBP_DATA_DIR/')

    return argument

def write_workflow_to_file(workflow_obj, filename):
    """
    This function writes the Broadband workflow contained in
    workflow_obj to the file specified in filename
    """
    install = install_cfg.InstallCfg.getInstance()

    workflowDoc = xml.getDOMImplementation().createDocument(None, "BBP_Run_Specification", None)
    root_element = workflowDoc.documentElement
    if workflow_obj.val_obj is not None:
        val_element = workflowDoc.createElement("Validation_Run")
        val_element.setAttribute("event", workflow_obj.val_obj.get_validation_name())
        val_element.setAttribute("version", install.VERSION)
        # Perform any substitutions we can
        workflow_obj.station_file = insert_variables(workflow_obj.station_file)
        val_element.setAttribute("input_station_file", workflow_obj.station_file)
        root_element.appendChild(val_element)
    else:
        val_element = workflowDoc.createElement("Scenario_Run")
        # Perform any substitutions we can
        workflow_obj.station_file = insert_variables(workflow_obj.station_file)
        val_element.setAttribute("input_station_file", workflow_obj.station_file)
        val_element.setAttribute("version", install.VERSION)
        root_element.appendChild(val_element)

    modules_element = workflowDoc.createElement("BBP_Modules")
    for mod in workflow_obj.workflow:
        mod_node = workflowDoc.createElement("BBP_Module")
        mod_name_node = workflowDoc.createElement("name")
        mod_name_node.appendChild(workflowDoc.createTextNode(mod.getName()))
        mod_node.appendChild(mod_name_node)
        files_node = workflowDoc.createElement("staged_files")
        for stage_file in mod.getStageFiles():
            file_node = workflowDoc.createElement("file")
            # Perform substitutions
            stage_file = insert_variables(stage_file)
            file_node.appendChild(workflowDoc.createTextNode(stage_file))
            files_node.appendChild(file_node)
        mod_node.appendChild(files_node)
        arguments_node = workflowDoc.createElement("arguments")
        for arg in mod.getArgs():
            #print arg.__class__
            if isinstance(arg, str):
                # Perform any substitutions we can
                arg = insert_variables(arg)
            arg_node = workflowDoc.createElement("argument")
            arg_node.setAttribute("type", arg.__class__.__name__)
            val = workflowDoc.createTextNode(str(arg))
            arg_node.appendChild(val)
            arguments_node.appendChild(arg_node)
        mod_node.appendChild(arguments_node)
        kw_args = mod.getKeywordArgs()
        if len(kw_args) > 0:
            kw_args_node = workflowDoc.createElement("keyword_arguments")
            for key in kw_args.keys():
                arg_node = workflowDoc.createElement("keyword_argument")
                arg_node.setAttribute("keyword", key)
                if isinstance(kw_args[key], validation_cfg.ValidationEvent):
                    arg_node.setAttribute("type", "validation_cfg.VE_EVENTS")
                    val = workflowDoc.createTextNode(kw_args[key].get_validation_name())
                    arg_node.appendChild(val)
                else:
                    arg_node.setAttribute("type", kw_args[key].__class__.__name__)
                    val = workflowDoc.createTextNode(str(kw_args[key]))
                    arg_node.appendChild(val)
                kw_args_node.appendChild(arg_node)
            mod_node.appendChild(kw_args_node)
        modules_element.appendChild(mod_node)
    root_element.appendChild(modules_element)
    out_fp = open(filename, "w")
    root_element.writexml(out_fp, indent="\t", addindent="\t", newl="\n")
    out_fp.close()

def write_workflow(workflow_obj, sim_id):
    """
    This function writes the Broadband workflow contained in
    workflow_obj to the XML directory using the current sim_id as
    filename.
    """
    install = install_cfg.InstallCfg.getInstance()
    write_workflow_to_file(workflow_obj, os.path.join(install.A_XML_DIR,
                                                      "%d.xml" % (sim_id)))

def parse_xml(xml_path):
    """
    This function parses a file with XML containing a Broadband
    workflow. It returns a workflow object.
    """
    workflow_list = []
    stat_file = ""
    version = ""
    val_obj = None
    dom_rep = xml.parse(xml_path)
    val_elements = dom_rep.getElementsByTagName("Validation_Run")
    if len(val_elements) > 0:
        val_element = val_elements[0]
        val_name = val_element.getAttribute("event")
        val_obj = getattr(validation_cfg.VE_EVENTS, val_name)
        version = val_element.getAttribute("version")

        stat_file = val_element.getAttribute("input_station_file")
        # Expand any variables that we find
        stat_file = expand_variables(stat_file)

    val_elements = dom_rep.getElementsByTagName("Scenario_Run")
    if len(val_elements) > 0:
        val_element = val_elements[0]
        version = val_element.getAttribute("version")
        stat_file = val_element.getAttribute("input_station_file")
        # Expand any variables that we find
        stat_file = expand_variables(stat_file)

    modules_element = dom_rep.getElementsByTagName("BBP_Modules")[0]
    # Print modules_element
    for mod in modules_element.childNodes:
        # BBP Module level
        if not mod.nodeType == mod.ELEMENT_NODE:
            continue
        # print mod
        module = Module()
        name_node = mod.getElementsByTagName("name")[0]
        name = name_node.childNodes[0].data.strip()
        # print name
        module.setName(name)
        files = mod.getElementsByTagName("staged_files")[0]
        for file_node in files.getElementsByTagName("file"):
            stage_file = file_node.childNodes[0].data.strip()
            # Expand any variables that we find
            stage_file = expand_variables(stage_file)

            module.addStageFile(stage_file)
        arguments = mod.getElementsByTagName("arguments")[0]
        for arg in arguments.getElementsByTagName("argument"):
            arg_type = arg.getAttribute("type")
            if not arg_type == "":
                parts = arg_type.split('.')
                if len(parts) > 1:
                    package = parts[0]
                    className = parts[1]
                    m = __import__(package)
                    objName = getattr(m, className)
                    arg_val = arg.childNodes[0].data.strip()
                    #print arg_val
                    arg_impl = getattr(objName, arg_val)
                else:
                    try:
                        m = __import__("__builtin__")
                    except ModuleNotFoundError:
                        m = __import__("builtins")
                    className = parts[0]
                    if className == "NoneType":
                        arg_impl = None
                    elif className == "dict":
                        arg_impl = ast.literal_eval(arg.childNodes[0].data.strip())
                    else:
                        objName = getattr(m, className)
                        # Handle the empty string case
                        if not arg.childNodes:
                            arg_val = ""
                        else:
                            arg_val = arg.childNodes[0].data.strip()
                        #print arg_val
                        arg_impl = objName(arg_val)
                if isinstance(arg_impl, str):
                    # Expand any variables that we have
                    arg_impl = expand_variables(arg_impl)
                module.addArg(arg_impl)
            else:
                arg_val = arg.childNodes[0].data.strip()
                module.addArg(arg_val)
        kw_args = mod.getElementsByTagName("keyword_arguments")
        if len(kw_args) > 0:
            for kw_arg in kw_args[0].getElementsByTagName("keyword_argument"):
                arg_keyword = str(kw_arg.getAttribute("keyword"))
                arg_type = kw_arg.getAttribute("type")
                if not arg_type == "":
                    parts = arg_type.split('.')
                    if len(parts) > 1:
                        package = parts[0]
                        className = parts[1]
                        m = __import__(package)
                        objName = getattr(m, className)
                        arg_val = kw_arg.childNodes[0].data.strip()
                        #print arg_val
                        arg_impl = getattr(objName, arg_val)
                    else:
                        try:
                            m = __import__("__builtin__")
                        except ModuleNotFoundError:
                            m = __import__("builtins")
                        className = parts[0]
                        if className == "NoneType":
                            arg_impl = None
                        else:
                            objName = getattr(m, className)
                            arg_val = kw_arg.childNodes[0].data.strip()
#                                                       print arg_val
                            if arg_type == "bool":
                                if str(arg_val).lower() == "false":
                                    arg_impl = False
                                else:
                                    arg_impl = True
                            else:
                                arg_impl = objName(arg_val)
#                                       print "kw: %s, val: %s, type: %s" % (arg_keyword, arg_impl, arg_type)
                    module.addKeywordArg(arg_keyword, arg_impl)

        workflow_list.append(module)
    workflow_obj = Workflow(workflow_list, val_obj, stat_file, version)
    return workflow_obj
