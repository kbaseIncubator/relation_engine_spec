.PHONY: test

test:
	docker-compose down
	docker-compose run spec sh /app/test/run_tests.sh
	docker-compose down