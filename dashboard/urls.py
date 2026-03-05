from django.urls import path
from django.views.generic import RedirectView
from . import views

urlpatterns = [
    path('', RedirectView.as_view(pattern_name='home', permanent=False)),
    path('home/', views.home, name='home'),
    path('strain/<str:strain_name>/', views.strain_view, name='strain'),
    path('mixed/', views.mixed_view, name='mixed'),
    path('mixed/pileup/', views.mixed_pileup_data, name='mixed_pileup_data'),

]
