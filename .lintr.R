linters <- linters_with_defaults(
  line_length_linter(256L),
  object_length_linter(64L),
  object_usage_linter = NULL,
  object_name_linter = NULL,
  brace_linter = NULL,
  spaces_left_parentheses_linter = NULL,
  commas_linter = NULL,
  paren_body_linter = NULL,
  commented_code_linter = NULL
)
