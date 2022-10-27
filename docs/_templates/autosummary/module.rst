{{ fullname | escape | underline}}

.. automodule:: {{ fullname }}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:
      :recursive:
   {% for item in functions %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:
      :recursive:
   {% for item in classes %}
      {{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
