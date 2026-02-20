from django.http import Http404
from django.shortcuts import render
from matplotlib import use

from .config import get_strain_config, get_strain_options
from .services import get_home_context, get_mixed_context, get_strain_context

use("agg")


def home(request):
    # Build page metrics/artifacts via service layer and keep view logic thin.
    context = get_home_context()
    context["strain_options"] = get_strain_options()
    return render(request, "home.html", context)


def strain_view(request, strain_name):
    # Only configured strains should resolve to a page.
    strain_config = get_strain_config(strain_name)
    if strain_config is None:
        raise Http404("Unknown strain")

    context = get_strain_context(strain_name, strain_config, request.get_host())
    context["strain_options"] = get_strain_options()
    return render(request, "strain.html", context)


def mixed_view(request):
    # Mixed page data flow is isolated in its own service module.
    context = get_mixed_context()
    context["strain_options"] = get_strain_options()
    return render(request, "mixed.html", context)
