from django.urls import path
from . import views

urlpatterns = [
    
    path('', views.analyze_sequence, name='analyze_sequence'),

   
    path('results/<int:analysis_id>/', views.results_view, name='results_view'),

    path('primer-designer/', views.primer_designer, name='primer_designer'),

    path('motif-finder/', views.motif_finder, name='motif_finder'),

    path("orf-finder/", views.orf_finder, name="orf_finder"),

    path("restriction-finder/", views.restriction_finder, name="restriction_finder"),

    
    path('sequence/list/', views.sequence_list, name='sequence_list'),   
    path('sequence/<int:pk>/', views.sequence_detail, name='sequence_detail'),  
    path('sequence/add/', views.sequence_create, name='sequence_create'),       
    path('sequence/<int:pk>/edit/', views.sequence_update, name='sequence_update'),  
    path('sequence/<int:pk>/delete/', views.sequence_delete, name='sequence_delete'),  
    path('save-sequence/', views.save_sequence, name='save_sequence'),
    path("saved-sequences/", views.saved_sequences, name="saved_sequences"),


]