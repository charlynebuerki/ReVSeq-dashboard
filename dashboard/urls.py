from django.urls import path
from . import views

urlpatterns = [
    path('home/', views.home, name='home'),
    path('strain/<str:strain_name>/', views.strain_view, name='strain'),
    path('mixed/', views.mixed_view, name='mixed'),

]