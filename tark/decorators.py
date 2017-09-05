from django.http import HttpResponse, HttpResponseRedirect, StreamingHttpResponse
import json
from tark.renderers import renderer, ALLOW_CONTENT_TYPE
import pprint
from tark.exceptions import FilterNotFound, AssemblyNotFound,\
    ReleaseNotFound
from __builtin__ import True

BOOL_PARAMETERS = ['filter_pk', 'expand_parent', 'expand', 'skip_sequence']

def render(function=None, default_content_type='application/json'):
    def decorator(view_func):
        def decorated(request, *args, **kwargs):
            
            try:
                response = view_func(request, *args, **kwargs)
            except (FilterNotFound, AssemblyNotFound, ReleaseNotFound) as e:
                response = renderer.render_error(str(e), **kwargs)
            except Exception as e:
                print "EXCEPTION: {}".format(type(e))
                print str(e)
                response = renderer.render_error("Unknown error", **kwargs)
            
            
            if isinstance(response, (HttpResponse, HttpResponseRedirect, StreamingHttpResponse)):
                return response
                        
            render_parameters = { 'expand': kwargs.get('expand', False),
                                 'expand_parent': kwargs.get('expand_parent', False),
                                 'filter_pk': True,
                                 'skip_sequence': kwargs.get('skip_sequence', False)}
    
            for media_type in request.accepted_types:
                if media_type.mimetype in ALLOW_CONTENT_TYPE:
                    render_parameters['content-type'] = media_type.mimetype

            if 'content-type' not in render_parameters:
                render_parameters['content-type'] = default_content_type

            return renderer.render(response, **render_parameters)
            
        decorated.__name__ = view_func.__name__
        decorated.__dict__ = view_func.__dict__
        decorated.__doc__ = view_func.__doc__

        return decorated
        
    if function is None:
        return decorator
    else:
        return decorator(function)
    
def parameter_parser(function=None, allow_methods=('GET', 'POST')):
    
    def decorator(view_func):
        def decorated(request, *args, **kwargs):
            
            if request.method == 'POST' and request.method in allow_methods:
                json_data = json.loads(request.body)
                kwargs.update(json_data)
            elif request.method == 'GET' and request.method in allow_methods:
                pass
            else:
                return HttpResponse(status=403)

            get_params = request.GET.dict()
            kwargs.update(get_params)

            for param in kwargs:
                if param in BOOL_PARAMETERS:
                    if kwargs[param] == '1' or kwargs[param] == 'True' or kwargs[param] == 'true':
                        kwargs[param] = True
                    elif kwargs[param] == '0' or kwargs[param] == 'False' or kwargs[param] == 'false':
                        kwargs[param] = False
                            
            response = view_func(request, *args, **kwargs)
            
            if isinstance(response, (HttpResponse, HttpResponseRedirect, StreamingHttpResponse)):
                return response   
            
        decorated.__name__ = view_func.__name__
        decorated.__dict__ = view_func.__dict__
        decorated.__doc__ = view_func.__doc__

        return decorated
        
    if function is None:
        return decorator
    else:
        return decorator(function)
    