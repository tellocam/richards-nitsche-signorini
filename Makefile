.PHONY: figures notebooks clean

figures:
	@for f in scripts/fig_test*.py; do uv run python "$$f"; done

notebooks:
	uv run jupyter lab notebooks/

clean:
	rm -f figures/*.pdf data/*.npz
