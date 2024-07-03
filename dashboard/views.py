from django.shortcuts import render
from django.http import HttpResponse
from django.template import loader
import pandas as pd
import geopandas as gpd
import pyproj
import os
from plotly.offline import plot
import plotly.express as px
import numpy as np
import plotly.graph_objects as go

color_dict = {
    "Adenovirus": "#65E0C8",
    "Bocavirus": "#D2D2E3",
    "Influenza A": "#6875D9",
    "Influenza B": "#E158D3",
    "Metapneumovirus": "royalblue",
    "Parainfluenza 2": "beige",
    "Parainfluenza 3": "#848877",
    "Parainfluenza 4": "#FDBF6F",
    "Polyomavirus": "#D9664D",
    "RSV - A/B": "thistle",
    "SARS-CoV-2" : "#E9D588",
    "Rhino- / Enterovirus": "#74457B",
    "coronavirus 229E": "pink",
    "coronavirus HKU1": "#BAE6B0",
    "coronavirus NL63": "#A968E1",
    "coronavirus OC43": "darkorange",
    "Mixed": "darkblue"
  }



def home(request):
    
    df = pd.read_csv('dashboard/static/example_data.csv')

    df['ent_date']= pd.to_datetime(df['ent_date'])
    df['week'] = df['ent_date'].map(lambda dt: dt.strftime('%Y/%W'))

    df['strain_name'][df['is_unique']==False] = 'Mixed'
    df['strains_PCR'] = pd.read_csv('dashboard/static/strains_PCR.tsv')['strains_PCR']

    df = df.drop_duplicates(subset='pseudonymized_id')
    

    grouped_df = df.groupby(['week','strain_name']).size().reset_index(name='count')

    fig = px.bar(grouped_df, x="week", y="count", color="strain_name",
            hover_data=['count'], barmode = 'stack',
            color_discrete_sequence=px.colors.qualitative.Light24,
            labels={
                     "week": "Week",
                     "count": "No. infections",
                     "strain_name": "Strain"
                 },
                title="Sequencing",
                color_discrete_map= color_dict,
                category_orders = {"strain_name": color_dict.keys()},
                template = "simple_white"
            )
    fig.update_xaxes(categoryorder='array', categoryarray= grouped_df['week'])
    fig.update_yaxes(showgrid=True)

    bar_seq = plot(fig, output_type="div")

    grouped_df = df.groupby(['week','strains_PCR']).size().reset_index(name='count')

    fig = px.bar(grouped_df, x="week", y="count", color="strains_PCR",
            hover_data=['count'], barmode = 'stack',
            color_discrete_sequence=px.colors.qualitative.Light24,
            labels={
                     "week": "Week",
                     "count": "No. infections",
                     "strains_PCR": "Strain"
                 },
            title="PCR",
            color_discrete_map= color_dict,
            category_orders = {"strains_PCR": color_dict.keys()},
            template = "simple_white"
            )
    fig.update_xaxes(categoryorder='array', categoryarray= grouped_df['week'])
    fig.update_yaxes(showgrid=True)

    bar_pcr = plot(fig, output_type="div")

    '''
    f="dashboard/static/swissBOUNDARIES3D_1_5_TLM_KANTONSGEBIET.shx"
    shapes = gpd.read_file(f, engine="pyogrio")

    shapes = shapes.to_crs('EPSG:4326')

    shapes['canton'] = ['GV', 'TG' ,'VS', 'AG', 'SZ', 'ZH' ,'OW', 'FR', 'GL', 'UR', 'NW', 'SO', 'AR', 'JU', 'GR', 'VD', 'LU', 'TI', 'ZG', 'BL', 'SG', 'SH', 'BE', 'BS', 'NE', 'AI']
    shapes=shapes[['geometry', 'canton', 'NAME']]

    shapes["geometry"] = (
        shapes.to_crs(shapes.estimate_utm_crs()).simplify(5000).to_crs(shapes.crs)
    )

    shapes.to_file('dashboard/static/swiss_cantons.geojson', driver='GeoJSON')
    

    geojson = gpd.read_file('dashboard/static/swiss_cantons.geojson')


    weeks = df['week'].unique()
    weeks.sort()

    cantons = np.array([[i]*len(weeks) for i in geojson.canton]).flatten()
    weeks = list(weeks) * len(geojson)

    data = pd.DataFrame({'canton': cantons, 'week': weeks})


    no_infections = df.groupby(['prescriber', 'week']).size().reset_index(name='count')
    no_infections.columns = ['canton', 'week', 'No. infections']
    
    data = data.merge(no_infections, how='left', on=['canton', 'week'])
    data['No. infections'] = data['No. infections'].fillna(0)

    data = geojson.merge(data, how='right', on= 'canton')


    fig = px.choropleth_mapbox(data,
                       geojson=data.geometry,
                       locations=data.index,
                       hover_name = 'NAME',
                       zoom = 5.5,
                       center = {"lat": 46.8, "lon": 8.5},
                       color='No. infections',
                       color_continuous_scale='viridis',
                       animation_frame='week',
                       mapbox_style = "carto-positron"
                       )
    texttrace = go.Scattermapbox(
        lat=geojson.geometry.centroid.y,
        lon=geojson.geometry.centroid.x,
        text=geojson["canton"],
        textfont={"color":"white","size":10},
        mode="text",
        name="canton",
        hoverinfo='skip'
    )

    fig.add_trace(texttrace)

    fig.write_html('dashboard/static/map.html')

    map_plt = plot(fig, output_type="div")
    '''
    
    return render(request, 'home.html', {'barplot_seq': bar_seq, 'barplot_pcr': bar_pcr }) 

def strain_view(request, strain_name):

    match strain_name:
        case 'Enterovirus':
            strain_name_long = 'Rhino- / Enterovirus'
        case 'RSV':
            strain_name_long = 'RSV - A/B'
        case _:
            strain_name_long= strain_name.replace('_' ,' ')
    
    #if not os.path.exists('dashboard/static/'+ strain_name):
    #   os.mkdir('dashboard/static/'+ strain_name)

    df = pd.read_csv('dashboard/static/example_data.csv')

    df['ent_date']= pd.to_datetime(df['ent_date'])
    df['week'] = df['ent_date'].map(lambda dt: dt.strftime('%Y/%W'))

    df = df[df['strain_name']==strain_name_long]

    if strain_name in ['Enterovirus', 'RSV']:
        if strain_name == 'Enterovirus':
            df['substrain_name'] = list(map({"Human rhinovirus A89": "Rhinovirus A", "Human enterovirus C109 isolate NICA08-4327": "Enterovirus C"}.get, df['substrain_name']))
        else:
            df['substrain_name'] = list(map({"Respiratory syncytial virus (type A)": "RSV A", "Human Respiratory syncytial virus 9320 (type B)": "RSV B"}.get, df['substrain_name']))

        grouped_df = df.groupby(['week', 'substrain_name']).size().reset_index(name='count')

        fig = px.bar(grouped_df, x="week", y="count", color="substrain_name",
            hover_data=['count'], barmode = 'stack',
            labels={
                     "week": "Week",
                     "count": "No. infections",
                     "substrain_name": "Substrain"
                 },
                color_discrete_map= color_dict,
                template = "simple_white"
            )
    

    else:

        grouped_df = df.groupby('week').size().reset_index(name='count')

        fig = px.bar(grouped_df, x="week", y="count", 
                labels={
                         "week": "Week",
                         "count": "No. infections"
                     },
                color_discrete_map= color_dict,
                template = "simple_white")
    

    fig.update_xaxes(categoryorder='array', categoryarray= grouped_df['week'])
    fig.update_yaxes(showgrid=True)

    bar_plt = plot(fig, output_type="div")

    geojson = gpd.read_file('dashboard/static/swiss_cantons.geojson')


    weeks = df['week'].unique()
    weeks.sort()

    cantons = np.array([[i]*len(weeks) for i in geojson.canton]).flatten()
    weeks = list(weeks) * len(geojson)

    data = pd.DataFrame({'canton': cantons, 'week': weeks})


    no_infections = df.groupby(['prescriber', 'week']).size().reset_index(name='count')
    no_infections.columns = ['canton', 'week', 'No. infections']
    
    data = data.merge(no_infections, how='left', on=['canton', 'week'])
    data['No. infections'] = data['No. infections'].fillna(0)

    data = geojson.merge(data, how='right', on= 'canton')

    if strain_name in ['Enterovirus', 'RSV']:
        substrain_df = pd.pivot_table(df, index=['prescriber','week'], columns='substrain_name', values='pseudonymized_id', aggfunc='count').reset_index()
        substrain_df.rename(columns={ 'prescriber': "canton" }, inplace = True)
        substrains = df['substrain_name'].unique()
        data = data.merge(substrain_df, how='left', on=['canton', 'week'])
        data = data.fillna(0)

        fig = px.choropleth_mapbox(data,
                       geojson=data.geometry,
                       locations=data.index,
                       hover_name = 'NAME',
                       hover_data = substrains,
                       zoom = 5.5,
                       center = {"lat": 46.8, "lon": 8.5},
                       color='No. infections',
                       color_continuous_scale='viridis',
                       animation_frame='week',
                       mapbox_style = "carto-positron"
                       )

    else:
        fig = px.choropleth_mapbox(data,
                           geojson=data.geometry,
                           locations=data.index,
                           hover_name = 'NAME',
                           zoom = 5.5,
                           center = {"lat": 46.8, "lon": 8.5},
                           color='No. infections',
                           color_continuous_scale='viridis',
                           animation_frame='week',
                           mapbox_style = "carto-positron"
                           )
    texttrace = go.Scattermapbox(
        lat=geojson.geometry.centroid.y,
        lon=geojson.geometry.centroid.x,
        text=geojson["canton"],
        textfont={"color":"white","size":10},
        mode="text",
        name="canton",
        hoverinfo='skip'
    )

    fig.add_trace(texttrace)
    map_plt = plot(fig, output_type="div")

    return render(request, "strain.html", {"strain_name": strain_name_long, 'map': map_plt ,'barplot': bar_plt}) 