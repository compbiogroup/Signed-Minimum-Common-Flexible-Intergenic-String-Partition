.PHONY: all test clean debug

all:
	stack build --test --no-run-tests --copy-bins --local-bin-path .

debug:
	stack build --fast --test --no-run-tests --copy-bins --local-bin-path .

test:
	stack build --fast --test

clean:
	stack clean
