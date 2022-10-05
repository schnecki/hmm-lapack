run-test:	update-test
	runhaskell Setup configure --user --enable-tests
	runhaskell Setup build
	runhaskell Setup haddock
	./dist/build/hmm-test/hmm-test

update-test:
	doctest-extract-0.1 -i src/ -i private/ -o test/ --library-main=Main --import-tested $$(cat test-module.list)
