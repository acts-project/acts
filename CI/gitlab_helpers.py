#! /usr/bin/env python
import requests

def test_response(response):
    response.raise_for_status()

def get_project_id(url,project_name,token):
    r = requests.get(url + "projects",headers={"PRIVATE-TOKEN":token})
    test_response(r)
    proj_list = r.json()
    matches = [proj_dict for proj_dict in proj_list if proj_dict["path_with_namespace"] == project_name]
    if len(matches) > 1:
        print "more than one project found matching the given project path"
        print "going to use the first one"
    return matches[0]["id"]

def get_merge_request_id(url,project_name,token,project_id,mr_iid):
    r = requests.get(url + "projects/{0}/merge_requests".format(project_id),headers={"PRIVATE-TOKEN":token})
    test_response(r)
    mr_list = r.json()
    matches = [mr_dict for mr_dict in mr_list if mr_dict["iid"] == mr_iid]
    if len(matches) > 1:
        print "more than one merge request found matching the given project path and merge request IID"
        print "going to use the first one"
    return matches[0]["id"]
    
