# Documentation Markdown Cheatsheet

Below there are some snippets for the Myst Markdown parser. You can find a full reference of the parser [here](https://myst-parser.readthedocs.io/en/latest/index.html).

## Pull in API documentation

* Code: 

````markdown
```{doxygenclass} Acts::Volume
---
members: center
---
```
````

* Result:

```{doxygenclass} Acts::Volume
---
members: center
---
```

(cheatsheetlabels)=
## Cross referencing and labels

* Setting a label:

```markdown
(cheatsheetlabels)=
```

* Referencing a label:

```markdown
Click [here](cheatsheetlabels) to come to the label.
```

Click [here](cheatsheetlabels) to come to the label.
