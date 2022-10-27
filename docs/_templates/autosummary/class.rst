{{ objname | escape | underline}}



.. autoclass:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ fullname }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}
   {% if methods %}
   .. rubric:: Methods

   .. autosummary::
   {% for item in methods %}
      {%- if item != '__init__' %}
      ~{{ fullname }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% for item in members %}
   {%- if item == 'Config' %}

Configuration object
####################

.. autoclass:: {{ fullname }}.Config
   :undoc-members:
   :members:
   {% endif %}
{%- endfor %}
