from django.http import Http404, JsonResponse
from django.shortcuts import render

from .config import get_strain_config, get_strain_options, is_mixed_page_enabled
from .services import get_home_context, get_mixed_context, get_mixed_pileup_context, get_strain_context

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

    context = get_strain_context(strain_name, strain_config)
    context["strain_options"] = get_strain_options()
    return render(request, "strain.html", context)


def mixed_view(request):
    if not is_mixed_page_enabled():
        raise Http404("Mixed page is disabled")
    # Mixed page data flow is isolated in its own service module.
    context = get_mixed_context(selected_sample=request.GET.get("sample"))
    context["strain_options"] = get_strain_options()
    return render(request, "mixed.html", context)


def mixed_pileup_data(request):
    if not is_mixed_page_enabled():
        raise Http404("Mixed page is disabled")
    payload = get_mixed_pileup_context(selected_sample=request.GET.get("sample"))
    return JsonResponse(payload)
